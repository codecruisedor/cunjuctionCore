#include "httplib.h"
#include <vector>
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include <iostream>
#include <chrono>
#include <libsgp4/SGP4.h>

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
    // Create an HTTP client pointing to the base domain
    httplib::Client cli("http://celestrak.org");

    // if (auto res = cli.Get("/NORAD/elements/gp.php?GROUP=planet&FORMAT=json"))
    // {
    //     const auto data = parseJson(res->body);
    //     // }
    // }
    // else
    // {
    //     std::cout << "GET request failed!" << std::endl;
    // }

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
            if (threeDone == 3) {
                threeDone = 0;
            }
            int j = i;
            while (j < res->body.size() && res->body[j] != '\n')
            {
                j++;
            }
            threeDone++;
           
            tleFirst = res->body.substr(i, j - i);
             if (threeDone == 1) {
                std::cout << "TLE First Line: " << tleFirst << std::endl;
            }else if (threeDone == 2) {
                tleSecond = res->body.substr(i, j - i);
                std::cout << "TLE Second Line: " << tleSecond << std::endl;
            }else if (threeDone == 3 ){
                tleThird = res->body.substr(i, j - i);
                std::cout << "TLE Third Line: " << tleThird << std::endl;
            }
            i = j + 1; // Move past the newline character
        }

        // std::cout << "Position (ECI): " << pos.ToString() << std::endl;

    }else {
        std::cout << "GET request failed!" << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    return 0;
}
