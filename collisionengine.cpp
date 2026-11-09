#include "httplib.h"
#include <vector>
#include <iostream>
#include <chrono>
#include <libsgp4/SGP4.h>
#include <libsgp4/Vector.h>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

std::set<std::pair<std::size_t, std::size_t>> unique_pairs;

struct Candidate {
    std::size_t first,second;
    std::size_t minIndex;
    double distance;
};
struct finalConjunctionpairs {
    std::size_t first,second;
    double missDistance;
    double relativeSpeed;
    std::string firstId, secondId;
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


inline double distSq(const libsgp4::Vector& a, const libsgp4::Vector& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}


inline double speed(const libsgp4::Vector& va, const libsgp4::Vector& vb) {
    const double dx = va.x - vb.x, dy = va.y - vb.y, dz = va.z - vb.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz); // km/s
}

// this packs the three cell coordinates into a single 64-bit integer, with an offset to allow negative coordinates
inline std::uint64_t cellKey(int cx, int cy, int cz) {
    constexpr int offset = 1 << 20;
    std::uint64_t cellkey = (std::uint64_t(cx + offset) << 42) | (std::uint64_t(cy + offset) << 21) | std::uint64_t(cz + offset);
    return cellkey;
}

inline std::uint64_t packPair(std::size_t i, std::size_t j) {  // requires i < j
    return (std::uint64_t(i) << 32) | std::uint64_t(j);
}

// Worst-case gate radius. Any pair's relative speed is at most 2 x the fastest
// object, so this never underestimates the gate used in screenPairs.
inline double maxGateKm(const std::vector<libsgp4::Vector>& velocities, const std::vector<char>& valid,
                        double stepSeconds, double thresholdKm, double marginKm) {
    double vmax = 0.0;
    for (std::size_t q = 0; q < velocities.size(); ++q) {
        if (!valid[q]) continue;
        const auto& v = velocities[q];
        vmax = std::max(vmax, std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z));
    }
    const double s = 2.0 * vmax * stepSeconds * 0.5;
    return std::sqrt(thresholdKm * thresholdKm + s * s) + marginKm;
}

// Every pair (i<j) that is within radiusKm at ANY grid sample.
std::unordered_set<std::uint64_t> candidatePairs(const std::vector<libsgp4::Vector>& positions,
                                                 const std::vector<char>& valid, std::size_t n,
                                                 std::size_t steps, double radiusKm) {
    std::unordered_set<std::uint64_t> pairs;
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> grid;
    std::vector<int> cx(n), cy(n), cz(n);
    const double r2 = radiusKm * radiusKm;

    for (std::size_t k = 0; k < steps; ++k) {
        grid.clear();
        for (std::size_t i = 0; i < n; ++i) {
            if (!valid[i * steps + k]) continue;

            const auto& p = positions[i * steps + k];

            cx[i] = static_cast<int>(std::floor(p.x / radiusKm));
            cy[i] = static_cast<int>(std::floor(p.y / radiusKm));
            cz[i] = static_cast<int>(std::floor(p.z / radiusKm));

            grid[cellKey(cx[i], cy[i], cz[i])].push_back(i);
        }
        for (std::size_t i = 0; i < n; ++i) {

            if (!valid[i * steps + k]) continue;

            const auto& pi = positions[i * steps + k];

                for (int dx = -1; dx <= 1; ++dx)
                for (int dy = -1; dy <= 1; ++dy)
                
                    for (int dz = -1; dz <= 1; ++dz) {
                        auto it = grid.find(cellKey(cx[i] + dx, cy[i] + dy, cz[i] + dz));
                        if (it == grid.end()) continue;
                        for (std::size_t j : it->second) {
                            if (j <= i) continue;                                  // each pair once
                            if (distSq(pi, positions[j * steps + k]) > r2) continue;
                            pairs.insert(packPair(i, j));
                        }
                    }
        }
    }
    return pairs;
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
                const double gate = std::sqrt(thr2 + ( halfTravel * halfTravel ) ) + gateMarginKm;

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


void screenOnePair(std::size_t i, std::size_t j,
                   const std::vector<libsgp4::Vector>& positions,
                   const std::vector<libsgp4::Vector>& velocities,
                   const std::vector<char>& valid, std::size_t steps,
                   double stepSeconds, double thr2, double gateMarginKm,
                   std::vector<double>& d2, std::vector<Candidate>& out)
{
    const libsgp4::Vector* ri = &positions[i * steps];
    const libsgp4::Vector* rj = &positions[j * steps];
    const libsgp4::Vector* vi = &velocities[i * steps];
    const libsgp4::Vector* vj = &velocities[j * steps];
    const char* oi = &valid[i * steps];
    const char* oj = &valid[j * steps];

    for (std::size_t k = 0; k < steps; ++k)
        d2[k] = (oi[k] && oj[k]) ? distSq(ri[k], rj[k]) : std::numeric_limits<double>::infinity();

    for (std::size_t k = 1; k + 1 < steps; ++k) {
        if (!(d2[k] <= d2[k - 1] && d2[k] < d2[k + 1])) continue;
        const double half = speed(vi[k], vj[k]) * stepSeconds * 0.5;
        const double gate = std::sqrt(thr2 + half * half) + gateMarginKm;
        if (d2[k] > gate * gate) continue;
        out.push_back({i, j, k, std::sqrt(d2[k])});
    }
}

std::optional<Conjunction> refineCandidate(
    const Candidate& c,
    const std::vector<libsgp4::SGP4>& propagators,
    const libsgp4::DateTime& start,
    double stepSeconds,
    double thresholdKm,
    double toleranceSeconds = 1e-3 ) {

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
    const double thresholdKm = 25.0; // Define the threshold for conjunctions
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<libsgp4::SGP4> propagators;
    std::vector<TLEData> tleDataList;
    std::unordered_map<std::string, std::int16_t> tleNameToIndex;
    std::vector<finalConjunctionpairs> highRiskCandidates;
    std::vector<finalConjunctionpairs> midRiskCandidates;
    std::vector<finalConjunctionpairs> lowRiskCandidates;

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
            //std::cout << "Satellite Name: " << t.name << std::endl;
            //std::cout << "Line 1: " << t.line1 << std::endl;
            //std::cout << "Line 2: " << t.line2 << std::endl;
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
        const double thr2 = thresholdKm * thresholdKm;
       // std::vector<double> d2(steps);
        const double radius = maxGateKm(velocities, valid, stepSeconds, thresholdKm, 10.0);
       // auto pairs = candidatePairs(positions, valid, n, steps, radius);
        auto startGrid = std::chrono::high_resolution_clock::now();
        auto pairs = candidatePairs(positions, valid, n, steps, radius);
        auto endGrid = std::chrono::high_resolution_clock::now();
        auto gridDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endGrid - startGrid);
        std::cout << "Grid screening time: " << gridDuration.count() << " milliseconds" << std::endl;
        std::cout << "Grid sample size: " << pairs.size() << std::endl;
        
        // std::cout << "Grid screening found " << filtered.size() << " candidate conjunctions" << std::endl;
        std::cout << "Propagation and screening completed." << std::endl;

        std::vector<double> d2(steps);
        std::vector<Candidate> filtered;

        for (std::uint64_t key : pairs) {
            const std::size_t i = key >> 32;               // top 32 bits
            const std::size_t j = key & 0xFFFFFFFFu;       // bottom 32 bits
            screenOnePair(i, j, positions, velocities, valid, steps,
                        stepSeconds, thr2, 10.0, d2, filtered);
                        
        }

        
        auto conjunctions = refineAll(filtered, propagators, start, stepSeconds, thresholdKm);
        std::cout << conjunctions.size() << " conjunctions found" << std::endl;
        for(const auto& conj : conjunctions) {
            unique_pairs.insert({conj.first, conj.second});
        }
        std::cout << "Total unique conjunctions found: " << unique_pairs.size() << " Out of total " << filtered.size() << std::endl;
        for(const auto& conj: conjunctions) {
            std::cout << "Conjunction between objects " << conj.first << " and " << conj.second << std::endl
            << "Miss distance: " << conj.missDistance << " km" << std::endl;
            if (conj.missDistance < 5.0) {
                highRiskCandidates.emplace_back(conj.first, conj.second, conj.missDistance, conj.relativeSpeed, tleDataList[conj.first].name, tleDataList[conj.second].name);
               // std::cout << "High risk conjunction!" << std::endl;
            } else if (conj.missDistance > 5.0 && conj.missDistance < 15.0) {
                midRiskCandidates.emplace_back(conj.first, conj.second, conj.missDistance, conj.relativeSpeed, tleDataList[conj.first].name, tleDataList[conj.second].name);
              //  std::cout << "Medium risk conjunction!" << std::endl;
            } else {
                lowRiskCandidates.emplace_back(conj.first, conj.second, conj.missDistance, conj.relativeSpeed, tleDataList[conj.first].name, tleDataList[conj.second].name);
               // std::cout << "Low risk conjunction." << std::endl;
            }
            //std::cout << "Conjunction between objects " << conj.first << " and " << conj.second << std::endl;
        }
        
        
    }
    else
    {
        std::cout << "GET request failed!" << std::endl;
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "High risk conjunctions: " << highRiskCandidates.size() << std::endl;
    std::cout << "Medium risk conjunctions: " << midRiskCandidates.size() << std::endl;
    std::cout << "Low risk conjunctions: " << lowRiskCandidates.size() << std::endl;
    std::cout << "Total unique conjunctions: " << unique_pairs.size() << std::endl;


    return 0;
}
