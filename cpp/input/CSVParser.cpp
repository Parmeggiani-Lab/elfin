#include "CSVParser.hpp"

#include "util.h"

namespace LokiAstari
{
std::string const& CSVRow::operator[](std::size_t index) const
{
	return m_data[index];
}
std::size_t CSVRow::size() const
{
	return m_data.size();
}
void CSVRow::readNextRow(std::istream& str)
{
	line = "";
	std::getline(str, line);

	std::stringstream   lineStream(line);
	std::string         cell;

	m_data.clear();
	while (std::getline(lineStream, cell, delimiter))
	{
		m_data.push_back(cell);
	}
	// This checks for a trailing comma with no data after it.
	if (!lineStream && cell.empty())
	{
		// If there was a trailing comma then add an empty element.
		m_data.push_back("");
	}
}
const std::string CSVRow::getLine() const
{
	return line;
}
int CSVRow::getNumParts() const
{
	return m_data.size();
}
const std::vector<std::string> & CSVRow::getParts() const
{
	return m_data;
}

std::istream& operator>>(std::istream& str, CSVRow& data)
{
	data.readNextRow(str);
	return str;
}

CSVIterator::CSVIterator(char _delimiter, std::istream& str) :
	m_str(str.good() ? & str : NULL),
	m_row(CSVRow(_delimiter))
{ ++(*this); }
CSVIterator::CSVIterator(std::istream& str) :
	m_str(str.good() ? & str : NULL)
{ ++(*this); }
CSVIterator::CSVIterator() :
	m_str(NULL)
{}

// Pre Increment
CSVIterator& CSVIterator::operator++()               {if (m_str) { if (!((*m_str) >> m_row)) {m_str = NULL;}} return *this;}
// Post increment
CSVIterator CSVIterator::operator++(int)             {CSVIterator    tmp(*this); ++(*this); return tmp;}
CSVRow const& CSVIterator::operator*()   const       {return m_row;}
CSVRow const* CSVIterator::operator->()  const       {return &m_row;}

bool CSVIterator::operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
bool CSVIterator::operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
}

namespace elfin
{

Points3f CSVParser::parseSpec(const std::string & filename)
{
	std::ifstream inputStream(filename);

	panic_if(!inputStream.is_open(),
	        "Could not open file: \"%s\"\n", filename.c_str());

	std::vector<float> data;
	for (LokiAstari::CSVIterator loop(' ', inputStream); loop != LokiAstari::CSVIterator(); ++loop)
	{
		const LokiAstari::CSVRow & row = (*loop);
		panic_if(
		    row.getNumParts() != 3, /* 3D points must be 3 values */
		    "Invalid row in spec file: \"%s\" - must have 3 components but there are %d\n",
		    row.getLine().c_str(), row.getNumParts());

		for (auto & part : row.getParts())
		{
			try {
				data.emplace_back(std::stof(part.c_str()));
			}
			catch (std::invalid_argument & e)
			{
				die("Failed to parse 3D float component at \"%s\" on function %s\n",
				    part.c_str(), e.what());
			}
		}
	}

	Points3f spec;

	for (std::vector<float>::iterator i = data.begin(); i != data.end(); i += 3)
		spec.emplace_back(i, i + 3);

	return spec;
}

} // namespace Elfin
