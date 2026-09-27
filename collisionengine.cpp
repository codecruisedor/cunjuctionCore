#include "httplib.h"
#include <vector>
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include <iostream>
#include <chrono>
#include <libsgp4/Tle.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Eci.h>

rapidjson::GenericArray<false, rapidjson::Value> parseJson(const std::string &jsonStr) {
    rapidjson::Document doc;
    doc.Parse(jsonStr.c_str());

    if (doc.HasParseError()) {
        std::cout << "Error parsing JSON" << std::endl;
        //return rapidjson::GenericArray<false, rapidjson::Value>(); // Return an empty array on parse error
    }

    //std::cout << "Number of elements in JSON array: " << doc.GetArray().Size() << std::endl;
    const auto &data = doc.GetArray();
    return data;
    

}

//void parseJson(const std::string &jsonStr);

int main() {
    auto start = std::chrono::high_resolution_clock::now();
      // Create an HTTP client pointing to the base domain
    httplib::Client cli("http://celestrak.org");

    if (auto res = cli.Get("/NORAD/elements/gp.php?GROUP=planet&FORMAT=json")) {
        // std::cout << "GET request successful!" << std::endl;
        // std::cout << "Response status: " << res->status << std::endl;
        //std::cout << "Response body: " << res->body << std::endl;

       const auto data = parseJson(res->body);
    //     for (const auto &item : data) {
    //         std::cout << "OBJECT_NAME: " << item["OBJECT_NAME"].GetString() << std::endl;
    //         // std::cout << "OBJECT_ID: " << item["OBJECT_ID"].GetString() << std::endl;
    //         // std::cout << "EPOCH: " << item["EPOCH"].GetString() << std::endl;
    //         // std::cout << "INCLINATION: " << item["INCLINATION"].GetDouble() << std::endl;
    //         // std::cout << "RA_OF_ASC_NODE: " << item["RA_OF_ASC_NODE"].GetDouble() << std::endl;
    //         // std::cout << "ECCENTRICITY: " << item["ECCENTRICITY"].GetDouble() << std::endl;
    //         // std::cout << "ARG_OF_PERICENTER: " << item["ARG_OF_PERICENTER"].GetDouble() << std::endl;
    //         // std::cout << "MEAN_ANOMALY: " << item["MEAN_ANOMALY"].GetDouble() << std::endl;
    //         // std::cout << "MEAN_MOTION: " << item["MEAN_MOTION"].GetDouble() << std::endl;
    //         // std::cout << "REV_AT_EPOCH: " << item["REV_AT_EPOCH"].GetInt() << std::endl;
    //         std::cout << "-----------------------------" << std::endl;
    // }
    } else {
        std::cout << "GET request failed!" << std::endl;
    }

    //actual TLE fetching
     if (auto res = cli.Get("/NORAD/elements/gp.php?GROUP=planet&FORMAT=TLE")) {
        std::cout << "GET request successful!" << std::endl;
        std::cout << "Response status: " << res->status << std::endl;
        std::cout << "Response body: " << res->body << std::endl;


    } else {
        std::cout << "GET request failed!" << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    return 0;
}






