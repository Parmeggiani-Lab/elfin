#ifndef _JSONPARSER_HPP_
#define _JSONPARSER_HPP_

#include "json.hpp"

#include "SpecParser.hpp"
#include "DBParser.hpp"
#include "../data/PairRelationship.hpp"

namespace elfin
{

// Credits to nolhmann from https://github.com/nlohmann/json
using JSON = nlohmann::json;

class JSONParser : public SpecParser, DBParser
{
public:
	JSONParser() {};
	virtual ~JSONParser() {};

	void parseDB(
	    const std::string & filename,
	    NameIdMap & nameMapOut,
	    IdNameMap & inmOut,
	    RelaMat & relMatOut,
	    RadiiList & radiiListOut);

	Points3f parseSpec(const std::string & filename);

	JSON parse(const std::string & filename);
	void inspect(const std::string & filename);
private:
};

} // namespace Elfin

#endif /* include guard */