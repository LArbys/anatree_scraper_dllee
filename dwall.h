#ifndef __DWALL_H__
#define __DWALL_H__

#include <vector>

namespace anascraper {
  
  float dwall( const std::vector<float>& pos, int& boundary_type );
  double dwall( const std::vector<double>& pos, int& boundary_type );

  float dspecificwall( const std::vector<float>& pos, const int boundary_type );
  double dspecificwall( const std::vector<double>& pos, const int boundary_type );
  
}

#endif
