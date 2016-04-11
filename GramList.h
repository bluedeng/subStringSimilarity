/*
 * GramList.h
 *
 *  Created on: Sep 14, 2010
 *      Author: xiaoliwang
 */
#ifndef _GREMLIST_H_
#define _GREMLIST_H_

#include <vector>
#include <utility>

using namespace std;

template <typename InvList = vector< pair<unsigned, int> > >
class CGramList {
 protected:
  vector< pair<unsigned, int> > invertedList;
  
 public:
  CGramList() {};
  InvList* getArray() { return &invertedList; }
  ~CGramList() {};
  void free() { delete this; }
  void clear() { }
};


class SPositionList {
protected:
	vector< pair<unsigned, unsigned> > positionList;

public:
	SPositionList() {};
	vector< pair<unsigned, unsigned> >* getArray() { return &positionList; }
	~SPositionList() {};
	void free() { delete this; }
	void clear() { }
};


#endif
