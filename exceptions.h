#pragma once
#include <exception>
#include <string>

class out_of_range_exception;
class generic_exception;
/** exceptions */
/**
 * @brief 下标越界异常
 *
 */
class out_of_range_exception : public std::exception
{
public:
	const char* what() {return "index out of range";}
};


/**
 * @brief 通用异常
 *
 */
class generic_exception : public std::exception
{
public:
	generic_exception(const char* w) throw()
	{
		msg = w;
	}
	virtual ~generic_exception() throw(){}
	const char* what() {return msg.c_str();}
private:
	std::string msg;
};
