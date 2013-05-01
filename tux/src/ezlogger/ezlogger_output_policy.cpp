#include "ezlogger_output_policy.hpp"

namespace axter
{

std::ostream& ezlogger_output_policy::get_log_stream()
{
  static const std::string FileName = EZLOGGER_OUTPUT_FILENAME;
#ifdef EZLOGGER_REPLACE_EXISTING_LOGFILE_
  static std::ofstream logfile(FileName.c_str(), std::ios_base::out);
#else
  static std::ofstream logfile(FileName.c_str(),  std::ios_base::out | std::ios_base::app);
#endif
  static bool logfile_is_open = logfile.is_open();
  if (logfile_is_open) return logfile;
  return std::cout;
}

}
