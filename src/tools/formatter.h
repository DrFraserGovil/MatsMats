#pragma once

#include <string>

/**
 * Convert all std::strings to const char* using constexpr if (C++17)
 */
template<typename T>
auto fmtconvert(T&& t) {
  if constexpr (std::is_same<std::remove_cv_t<std::remove_reference_t<T>>, std::string>::value) {
    return std::forward<T>(t).c_str();
  }
  else {
    return std::forward<T>(t);
  }
}

/**
 * printf like formatting for C++ with std::string
 * Original source: https://stackoverflow.com/a/26221725/11722
 */
template<typename ... Args>
std::string stringFormatInternal(const std::string& format, Args&& ... args)
{
  const auto size = snprintf(nullptr, 0, format.c_str(), std::forward<Args>(args) ...) + 1;
  if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
  std::unique_ptr<char[]> buf(new char[size]);
  snprintf(buf.get(), size, format.c_str(), args ...);
  return std::string(buf.get(), buf.get() + size - 1);
}

template<typename ... Args>
std::string stringFormat(std::string fmt, Args&& ... args) {
  return stringFormatInternal(fmt, fmtconvert(std::forward<Args>(args))...);
}


// Function to format directly into an existing std::string object.
// This uses the two-pass snprintf approach for robust sizing.
template<typename ... Args>
void format_to_string(std::string& out_str, const std::string& format, Args&& ... args) {

	const auto size = snprintf(nullptr, 0, format.c_str(), fmtconvert(std::forward<Args>(args))...) + 1;
	if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
	out_str.resize(size - 1);
	const int written_chars = snprintf(out_str.data(), size, format.c_str(), fmtconvert(std::forward<Args>(args))...);
	 if (written_chars < 0 || written_chars >= size) {
        // This indicates an error during formatting or a buffer overflow (shouldn't happen if required_size_with_null was correct)
        throw std::runtime_error("Error during string formatting: snprintf failed or truncated in second pass.");
    }
	 out_str.resize(written_chars);	
}