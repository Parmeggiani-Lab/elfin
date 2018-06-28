#include <fstream>

#include "JSONParser.hpp"
#include "util.h"

namespace elfin
{

void JSONParser::parseDB(
    const std::string & filename,
    NameIdMap & nameMapOut,
    IdNameMap & inmOut,
    RelaMat & relMatOut,
    RadiiList & radiiListOut)
{
	JSON j = this->parse(filename);
	if (relMatOut.size() > 0)
		wrn("JSONParser::parseDB(): argument relMatOut is not empty!");

	relMatOut.clear();
	JSON doublesData = j["doublesData"];
	const size_t dim = doublesData.size();
	relMatOut.resize(dim);

	// First we need the name-to-ID map
	unsigned int id = 0;
	for (JSON::iterator it = doublesData.begin();
	        it != doublesData.end();
	        ++it)
	{
		// JSON::iterator does not support offsetting,
		// but we need the key() which is only available
		// from JSON::iterator
		nameMapOut[it.key()] = id;
		inmOut[id] = it.key();

		RelaRow & row = relMatOut.at(id);
		row.resize(dim, NULL);

		dbg("Name-ID Pair: %s<->%d, neighbours: %d\n",
		    it.key().c_str(), id, (*it).size());
		id++;
	}

	// With the mame-to-ID map we can construct the
	// relationship matrix without using string names
	//
	// This is thinking ahead for the data to be operated
	// on by a more C-like kernel. In that case strings
	// strings are definitely bad for performance
	for (JSON::iterator it = doublesData.begin();
	        it != doublesData.end();
	        ++it) {
		id = nameMapOut[it.key()];
		RelaRow & row = relMatOut.at(id);

		// Add neighbouring nodes
		for (JSON::iterator innerIt = (*it).begin();
		        innerIt != (*it).end();
		        ++innerIt) {
			const unsigned int innerId = nameMapOut[innerIt.key()];
			prf("innerIt key: %s, id: %u, data: %s\n\n",
			    innerIt.key().c_str(),
			    innerId,
			    innerIt->dump().c_str());

			panic_if(
			    row.at(innerId) != NULL,
			    "Initial PairRelationship must be NULL - implementation error?\n");

			std::vector<float> combBv;
			for (JSON::iterator comBf = (*innerIt)["comB"].begin();
			        comBf != (*innerIt)["comB"].end();
			        ++comBf) {
				combBv.push_back((*comBf).get<float>());
			}
			std::vector<float> rotv;
			for (JSON::iterator rotRow = (*innerIt)["rot"].begin();
			        rotRow != (*innerIt)["rot"].end();
			        ++rotRow) {
				for (JSON::iterator rotF = (*rotRow).begin();
				        rotF != (*rotRow).end();
				        ++rotF) {
					rotv.push_back((*rotF).get<float>());
				}
			}
			std::vector<float> tranv;
			for (JSON::iterator tranRow = (*innerIt)["tran"].begin();
			        tranRow != (*innerIt)["tran"].end();
			        ++tranRow) {
				for (JSON::iterator tranF = (*tranRow).begin();
				        tranF != (*tranRow).end();
				        ++tranF) {
					tranv.push_back((*tranF).get<float>());
				}
			}
			row.at(innerId) = new PairRelationship(
			    combBv, rotv, tranv);
		}

		id++;
	}

	// Verify parsed pair relationship
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			PairRelationship * pr = relMatOut.at(i).at(j);
			if (pr != NULL)
			{
				prf("relMatOut[%2d][%2d] is:\n%s\n",
				    i, j, pr->toString().c_str());
			}
		}
	}

	// Parse radii for collision checking
	JSON singlesData = j["singlesData"];
	for (JSON::iterator it = singlesData.begin();
	        it != singlesData.end();
	        ++it)
	{
		JSON & radii = (*it)["radii"];
		radiiListOut.emplace_back(radii["avgAll"],
		                          radii["maxCA"],
		                          radii["maxHeavy"]);
	}
}

Points3f JSONParser::parseSpec(
    const std::string & filename)
{
	JSON j = this->parse(filename);

	std::vector<float> data;
	for (auto & point3d : j["coms"])
		for (auto & part3d : point3d)
			data.push_back(part3d.get<float>());

	panic_if(data.size() % 3 != 0,
	         "Input data not in the form of 3D tuple values!\n");

	Points3f spec;

	for (std::vector<float>::iterator i = data.begin(); i != data.end(); i += 3)
		spec.emplace_back(i, i + 3);

	return spec;
}

JSON JSONParser::parse(const std::string & filename)
{
	std::ifstream inputStream(filename);

	panic_if(!inputStream.is_open(),
	         "Could not open file: \"%s\"\n", filename.c_str());

	JSON j;
	inputStream >> j;

	return j;
}

void JSONParser::inspect(const std::string & filename)
{
	std::cout << this->parse(filename).dump() << std::endl;
}

} // namespace elfin