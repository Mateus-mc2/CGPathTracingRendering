#ifndef UTIL_EXCEPTION_H_
#define UTIL_EXCEPTION_H_

#include <string>

namespace util {

class UtilException : public std::exception {
  private:
    const std::string kErrorMsg;
  public:
    explicit UtilException(const std::string &error) : kErrorMsg(error) {}

    const char* what() const { return this->kErrorMsg.c_str(); }
};

}  // namespace util

#endif  // UTIL_EXCEPTION_H_