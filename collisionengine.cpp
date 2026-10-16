#include "httplib.h"
#include <vector>
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include <iostream>
#include <chrono>
#include <libsgp4/SGP4.h>
#include <libsgp4/Vector.h>
#include <cmath>
#include <cstddef>

void createEarth() {

}
struct Candidate {
    std::size_t first,second;
    std::size_t minIndex;
    double distance;
};
struct Conjunction {
    std::size_t first,second;
    libsgp4::DateTime timeOfClosestApproach;
    double missDistance;
    double relativeSpeed;
};
struct TLEData
{
    std::string name;
    std::string line1;
    std::string line2;
};

rapidjson::GenericArray<false, rapidjson::Value> parseJson(const std::string &jsonStr)
{
    rapidjson::Document doc;
    doc.Parse(jsonStr.c_str());

    if (doc.HasParseError())
    {
        std::cout << "Error parsing JSON" << std::endl;
        // return rapidjson::GenericArray<false, rapidjson::Value>(); // Return an empty array on parse error
    }

    // std::cout << "Number of elements in JSON array: " << doc.GetArray().Size() << std::endl;
    const auto &data = doc.GetArray();
    return data;
}

// void parseJson(const std::string &jsonStr);
inline double distSq(const libsgp4::Vector& a, const libsgp4::Vector& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

inline double speed(const libsgp4::Vector& va, const libsgp4::Vector& vb) {
    const double dx = va.x - vb.x, dy = va.y - vb.y, dz = va.z - vb.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz); // km/s
}

std::vector<Candidate> screenPairs(
    const std::vector<libsgp4::Vector>& positions,
    const std::vector<libsgp4::Vector>& velocities,
    const std::vector<char>& valid,
    std::size_t n,
    std::size_t steps,
    double stepSeconds,
    double thresholdKm,
    double gateMarginKm = 10.0)
{
    std::vector<Candidate> candidates;
    std::vector<double> d2(steps); // reused per pair, avoids reallocating

    const double thr2 = thresholdKm * thresholdKm;

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const libsgp4::Vector* ri = &positions[i * steps];
            const libsgp4::Vector* rj = &positions[j * steps];
            const char* oi = &valid[i * steps];
            const char* oj = &valid[j * steps];

            // build this pair's distance^2 curve from the cached grid
            for (std::size_t k = 0; k < steps; ++k) {
                d2[k] = (oi[k] && oj[k]) ? distSq(ri[k], rj[k]) : std::numeric_limits<double>::infinity();
            }

            // scan for sampled local minima
            for (std::size_t k = 1; k + 1 < steps; ++k) {
                if (!(d2[k] <= d2[k - 1] && d2[k] < d2[k + 1])) 
                    continue;

                // gate: is this dip worth refining?
                const libsgp4::Vector* vi = &velocities[i * steps];
                const libsgp4::Vector* vj = &velocities[j * steps];
                const double vrel = speed(vi[k], vj[k]);
                const double halfTravel = vrel * stepSeconds * 0.5;
                const double gate = std::sqrt(thr2 + halfTravel * halfTravel) + gateMarginKm;

                if (d2[k] > gate * gate) continue; // safe to discard

                candidates.push_back({i, j, k, std::sqrt(d2[k])});
            }
        }
    }
    return candidates;
}

inline bool stateAt(const libsgp4::SGP4& sgp4, const libsgp4::DateTime& t,
                    libsgp4::Vector& r, libsgp4::Vector& v) {
    try {
        const libsgp4::Eci eci = sgp4.FindPosition(t);
        r = eci.Position();
        v = eci.Velocity();
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

// Squared distance between two objects at `secondsFromStart`.
// Returns infinity if either propagation fails.
inline double pairDistSq(const libsgp4::SGP4& a, const libsgp4::SGP4& b,
                         const libsgp4::DateTime& start, double secondsFromStart) {
    libsgp4::Vector rX, vY, rZ, vW;
    const libsgp4::DateTime t = start.AddSeconds(secondsFromStart);
    if (!stateAt(a, t, rX, vY) || !stateAt(b, t, rZ, vW))
        return std::numeric_limits<double>::infinity();
    return distSq(rX, rZ);  // the helper from screenPairs
}

std::optional<Conjunction> refineCandidate(
    const Candidate& c,
    const std::vector<libsgp4::SGP4>& propagators,
    const libsgp4::DateTime& start,
    double stepSeconds,
    double thresholdKm,
    double toleranceSeconds = 1e-3)
{
    const libsgp4::SGP4& A = propagators[c.first];
    const libsgp4::SGP4& B = propagators[c.second];

    constexpr double invPhi = 0.6180339887498949;  // 1 / golden ratio

    // Bracket: the samples on either side of the sampled minimum
    double lo = static_cast<double>(c.minIndex - 1) * stepSeconds;
    double hi = static_cast<double>(c.minIndex + 1) * stepSeconds;

    double x1 = hi - invPhi * (hi - lo);
    double x2 = lo + invPhi * (hi - lo);
    
    double f1 = pairDistSq(A, B, start, x1);
    double f2 = pairDistSq(A, B, start, x2);

    while (hi - lo > toleranceSeconds) {
        if (f1 < f2) {              // minimum is in [lo, x2]
            hi = x2;
            x2 = x1;
            f2 = f1;       // reuse the old probe
            x1 = hi - invPhi * (hi - lo);
            f1 = pairDistSq(A, B, start, x1);
        } else {                    // minimum is in [x1, hi]
            lo = x1;
            x1 = x2;
            f1 = f2;       // reuse the old probe
            x2 = lo + invPhi * (hi - lo);
            f2 = pairDistSq(A, B, start, x2);
        }
    }

    const double sTca = 0.5 * (lo + hi);
    const libsgp4::DateTime tca = start.AddSeconds(sTca);

    libsgp4::Vector ra, va, rb, vb;
    if (!stateAt(A, tca, ra, va) || !stateAt(B, tca, rb, vb))
        return std::nullopt;

    const double miss = std::sqrt(distSq(ra, rb));
    if (miss > thresholdKm) return std::nullopt;  // the strict test

    return Conjunction{c.first, c.second, tca, miss, speed(va, vb)};
}

std::vector<Conjunction> refineAll(
    const std::vector<Candidate>& candidates,
    const std::vector<libsgp4::SGP4>& propagators,
    const libsgp4::DateTime& start,
    double stepSeconds,
    double thresholdKm)
{
    std::vector<Conjunction> out;
    for (const auto& c : candidates) {
        if (auto conj = refineCandidate(c, propagators, start, stepSeconds, thresholdKm))
            out.push_back(*conj);
    }
    std::sort(out.begin(), out.end(),
              [](const Conjunction& a, const Conjunction& b) { return a.missDistance < b.missDistance; });
    return out;
}

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<libsgp4::SGP4> propagators;
    std::vector<TLEData> tleDataList;
    std::unordered_map<std::string, std::int16_t> tleNameToIndex;

    // Create an HTTP client pointing to the base domain
    httplib::Client cli("http://celestrak.org");

    // actual TLE fetching
    if (auto res = cli.Get("/NORAD/elements/gp.php?GROUP=satnogs&FORMAT=TLE"))
    {
        std::string tleFirst = "";
        std::string tleSecond = "";
        std::string tleThird = "";
        int threeDone = 0;
        int router = 0;
        std::cout << "GET request successful!" << std::endl;
        std::cout << "Response status: " << res->status << std::endl;
        // std::cout << "Response body: " << res->body << std::endl;
        // for (int u = 0; u < res->body.size(); u++)
        // {
        //     // if (res->body[u] == '\n')
        //     // {
        //     //     std::cout << "Newline found at index: " << u << std::endl;
        //     // }
        // }
        for (int i = 0; i < res->body.size();) {
            if (threeDone == 3){
                 threeDone = 0;
            }
            int j = i;
            while (j < res->body.size() && res->body[j] != '\n'){
                j++;
            }   
            threeDone++;
            if (threeDone == 1){
                tleFirst  =  res->body.substr(i, j - i);
            }
            else if (threeDone == 2){
                tleSecond =  res->body.substr(i, j - i);
            }    
            else if (threeDone == 3)
            {
                tleThird = res->body.substr(i, j - i);
                // Create a TLEData object and add it to the list
                TLEData tleData;
                if (tleNameToIndex.find(tleFirst) != tleNameToIndex.end()) {
                    std::cout << "Duplicate TLE name found: " << tleFirst << std::endl;
                   // continue; // Skip this TLE and move to the next one
                } else {
                    tleNameToIndex[tleFirst] = static_cast<std::int16_t>(tleDataList.size());
                    tleData.name = tleFirst;
                    tleData.line1 = tleSecond.substr(0,69);
                    tleData.line2 = tleThird.substr(0,69);
                    tleDataList.push_back(tleData);
                }
            }
            i = j + 1; // Move past the newline character
        }

        // iterate through the TLE data list and print the TLEs
        auto startProp = std::chrono::high_resolution_clock::now();
        std::cout << "SIZE : " << tleDataList.size();
        
        for (const auto &t : tleDataList) {
            std::cout << "Satellite Name: " << t.name << std::endl;
            std::cout << "Line 1: " << t.line1 << std::endl;
            std::cout << "Line 2: " << t.line2 << std::endl;
            libsgp4::Tle tle(t.name, t.line1, t.line2); 
            libsgp4::SGP4 sgp4(tle);
            propagators.emplace_back(tle);
        }

        // shared window: same start time and step for every object
        libsgp4::DateTime start = libsgp4::DateTime::Now(true); // true = UTC
        const double stepSeconds = 5.0;
        const double windowHours = 24.0;
        const std::size_t steps = static_cast<std::size_t>(windowHours * 3600.0 / stepSeconds) + 1;

        const std::size_t n = propagators.size();
        std::cout << "Propagating " << n << " objects for " << steps << " steps each" << std::endl;
        int totalPropgationSteps = n * steps;
        std::vector<libsgp4::Vector> positions(totalPropgationSteps);
        std::vector<libsgp4::Vector> velocities(totalPropgationSteps);
        std::vector<char> valid(totalPropgationSteps, 0);

        
       for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 0; k < steps; ++k) {

            libsgp4::DateTime t = start.AddSeconds(static_cast<double>(k) * stepSeconds);
            try {
                libsgp4::Eci eci = propagators[i].FindPosition(t);
                positions[i * steps + k]  = eci.Position();
                velocities[i * steps + k] = eci.Velocity();
                valid[i * steps + k] = 1;
                //std::cout << "Satellite " << i << ", Step " << k << ": Position (ECI): " << positions[i * steps + k].ToString() << std::endl;
                }catch (const std::exception&) {
                    valid[i * steps + k] = 0; // decayed orbit, bad propagation, etc.
                }
            }
        }
        auto candidates = screenPairs(positions, velocities, valid, n, steps, stepSeconds, 10.0);
        std::cout << candidates.size() << " candidate conjunctions found" << std::endl;

        auto endProp = std::chrono::high_resolution_clock::now();
        auto propDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endProp - startProp);
        std::cout << "Propagation time: " << propDuration.count() << " milliseconds" << std::endl;
        std::cout << "Propagation position size : " << totalPropgationSteps << std::endl;
        const double thresholdKm = 25.0; // Define the threshold for conjunctions

        //call refineAll to get the conjunctions
        auto conjunctions = refineAll(candidates, propagators, start, stepSeconds, thresholdKm);

        std::cout << conjunctions.size() << " conjunctions found" << std::endl;
        for (const auto& conj : conjunctions) {
            std::cout << "Conjunction between objects " << conj.first << " and " << conj.second << " at time " << conj.timeOfClosestApproach.ToString() << " with miss distance " << conj.missDistance << " km and relative speed " << conj.relativeSpeed << " km/s" << std::endl;
        }
        // std::cout << "Position (ECI): " << pos.ToString() << std::endl;
    }
    else
    {
        std::cout << "GET request failed!" << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;

   
   
    return 0;
}
