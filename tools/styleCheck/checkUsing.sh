# check of no using's in .hpp files
grep -rn using ../../scai | grep hpp | grep -v "\*" | grep -v "//" | grep -v "::"

#only known exception
#scai/logging.hpp:199:#define SCAI_LOG_USING(alogger) using alogger;
