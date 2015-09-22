# - Try to find Catch
# Once done this will define
#  FOUND_CATCH - If Catch headers were found
#  Catch_INCLUDE_DIRS - The Catch include directory

find_path(CATCH_INCLUDE_DIR catch.hpp)

if (CATCH_INCLUDE_DIR)
  set(FOUND_CATCH TRUE)
  message(STATUS "Found header for Catch")
endif (CATCH_INCLUDE_DIR)

mark_as_advanced(CATCH_INCLUDE_DIR) 
