#ifndef _CSVPARSER_HPP_
#define _CSVPARSER_HPP_

#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>
#include <fstream>
#include <string>

#include "SpecParser.hpp"

/*
 * A Nice CSV Parser Iterator by Loki Astari from Stack Overflow
 * Credits to http://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
 * Some modifications have been made
 *
 * Note: not that I did not know how to write a CSV parser, but it's
 * nice to learn from what's generally received as a 'nice and clean'
 * implementation, especially in a style far from how C does it
 *
 */

namespace LokiAstari
{

class CSVRow
{
public:
	CSVRow() : delimiter(',') {}
	CSVRow(char _delimiter) : delimiter(_delimiter) {}
	std::string const& operator[](std::size_t index) const;
	std::size_t size() const;
	void readNextRow(std::istream& str);
	const std::string getLine() const;
	int getNumParts() const;
	const std::vector<std::string> & getParts() const;
private:
	const char delimiter;
	std::vector<std::string> m_data;
	std::string line;
};

class CSVIterator
{
public:
	typedef std::input_iterator_tag     iterator_category;
	typedef CSVRow                      value_type;
	typedef std::size_t                 difference_type;
	typedef CSVRow*                     pointer;
	typedef CSVRow&                     reference;

	CSVIterator(char _delimiter, std::istream& str);
	CSVIterator(std::istream& str);
	CSVIterator();

	CSVIterator& operator++();
	CSVIterator operator++(int);
	CSVRow const& operator*() const;
	CSVRow const* operator->() const;

	bool operator==(CSVIterator const& rhs);
	bool operator!=(CSVIterator const& rhs);
private:
	std::istream*       m_str;
	CSVRow              m_row;
};

} // namespace LokiAstari

namespace elfin
{

class CSVParser : public SpecParser
{
public:
	CSVParser() {};
	virtual ~CSVParser() {};

	Points3f parseSpec(const std::string & filename);

private:
};

} // namespace CSVParser

#endif /* include guard */