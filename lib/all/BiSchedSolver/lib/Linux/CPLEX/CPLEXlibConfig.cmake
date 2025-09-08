# CPLEX Linux edition

set(CPLEX_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib/Linux/CPLEX/cplex/include/;${CMAKE_CURRENT_SOURCE_DIR}/lib/Linux/CPLEX/concert/include/;${CMAKE_CURRENT_SOURCE_DIR}/lib/Linux/CPLEX/cpoptimizer/include/")

string(STRIP "${CPLEX_INCLUDE_DIRS}" CPLEX_INCLUDE_DIRS)

# use cplex lib /!\ the order is important
set(CPLEX_LIBRARIES "${CMAKE_CURRENT_SOURCE_DIR}/lib/Linux/CPLEX/cplex/lib/x86-64_linux/static_pic/libilocplex.a;${CMAKE_CURRENT_SOURCE_DIR}/lib/Linux/CPLEX/cplex/lib/x86-64_linux/static_pic/libcplex.a;${CMAKE_CURRENT_SOURCE_DIR}/lib/Linux/CPLEX/concert/lib/x86-64_linux/static_pic/libconcert.a;${CMAKE_CURRENT_SOURCE_DIR}/lib/Linux/CPLEX/cpoptimizer/lib/x86-64_linux/static_pic/libcp.a")

string(STRIP "${CPLEX_LIBRARIES}" CPLEX_LIBRARIES)
