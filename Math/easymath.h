#pragma once

#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <utility>
#include <set>
#include <list>

namespace easymath{
	class XY{
	public:
		XY(double x, double y):x(x),y(y){};
		XY(){};
		~XY(){};
		double x,y;
		XY operator-(const XY &other){
			return XY(x-other.x, y-other.y);
		}
		friend bool operator<(const XY &lhs, const XY &rhs){
			if (lhs.x!=rhs.x) return lhs.x<rhs.x;
			return lhs.y<rhs.y; // if same, sort by y values
		}
		friend bool operator==(const XY &lhs, const XY &rhs){
			return lhs.x==rhs.x && lhs.y==rhs.y;
		}
	};

}

