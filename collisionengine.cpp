#include "httplib.h"
#include <vector>
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include <iostream>
#include <chrono>
#include <libsgp4/SGP4.h>

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

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<libsgp4::SGP4> propagators;
    std::vector<TLEData> tleDataList;

    // Create an HTTP client pointing to the base domain
    httplib::Client cli("http://celestrak.org");

    // actual TLE fetching
    if (auto res = cli.Get("/NORAD/elements/gp.php?GROUP=planet&FORMAT=TLE"))
    {
        std::string tleFirst = "";
        std::string tleSecond = "";
        std::string tleThird = "";
        int threeDone = 0;
        int router = 0;
        std::cout << "GET request successful!" << std::endl;
        std::cout << "Response status: " << res->status << std::endl;
        // std::cout << "Response body: " << res->body << std::endl;
        for (int u = 0; u < res->body.size(); u++)
        {
            if (res->body[u] == '\n')
            {
                std::cout << "Newline found at index: " << u << std::endl;
            }
        }
        for (int i = 0; i < res->body.size();)
        {
            if (threeDone == 3)
            {
                threeDone = 0;
            }
            int j = i;
            while (j < res->body.size() && res->body[j] != '\n')
            {
                j++;
            }
            threeDone++;
            if (threeDone == 1)
            {
                tleFirst = res->body.substr(i, j - i);
               // std::cout << "TLE First Line: " << tleFirst << std::endl;
            }
            else if (threeDone == 2)
            {
                tleSecond = res->body.substr(i, j - i);
               // std::cout << "TLE Second Line: " << tleSecond << std::endl;
            }
            else if (threeDone == 3)
            {
                tleThird = res->body.substr(i, j - i);
               // std::cout << "TLE Third Line: " << tleThird << std::endl;
                // Create a TLEData object and add it to the list
                TLEData tleData;
                tleData.name = tleFirst; // You can set this to a meaningful name if available
                tleData.line1 = tleSecond.substr(0,69);
                tleData.line2 = tleThird.substr(0,69);
                tleDataList.push_back(tleData);
            }
            i = j + 1; // Move past the newline character
        }

        // iterate through the TLE data list and print the TLEs
        auto startProp = std::chrono::high_resolution_clock::now();
        std::cout << "SIZE : " << tleDataList.size();
        
        for (const auto &t : tleDataList) {
            libsgp4::Tle tle(t.name, t.line1, t.line2); 
            libsgp4::SGP4 sgp4(tle);
            propagators.emplace_back(tle);
        }

        // shared window: same start time and step for every object
        libsgp4::DateTime start = libsgp4::DateTime::Now(true); // true = UTC
        const double stepSeconds = 30.0;
        const double windowHours = 24.0;
        const std::size_t steps = static_cast<std::size_t>(windowHours * 3600.0 / stepSeconds) + 1;

        const std::size_t n = propagators.size();

        std::vector<libsgp4::Vector> positions(n * steps);
        std::vector<libsgp4::Vector> velocities(n * steps);
        std::vector<char> valid(n * steps, 0);

        
       
       for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 0; k < steps; ++k) {

            libsgp4::DateTime t = start.AddSeconds(static_cast<double>(k) * stepSeconds);
            try {
                libsgp4::Eci eci = propagators[i].FindPosition(t);
                positions[i * steps + k]  = eci.Position();
                velocities[i * steps + k] = eci.Velocity();
                valid[i * steps + k] = 1;
                std::cout << "Satellite " << i << ", Step " << k << ": Position (ECI): " << positions[i * steps + k].ToString() << std::endl;
                }catch (const std::exception&) {
                    valid[i * steps + k] = 0; // decayed orbit, bad propagation, etc.
                }
            }
        }
        auto endProp = std::chrono::high_resolution_clock::now();
        auto propDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endProp - startProp);
        std::cout << "Propagation time: " << propDuration.count() << " milliseconds" << std::endl;
         std::cout << "Propagation position size : " << n * steps << std::endl;

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
