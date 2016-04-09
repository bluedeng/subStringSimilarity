/*
 * SeqDB.h
 *
 *  Created on: Sep 14, 2010
 *      Author: xiaoliwang
 */

#ifndef _SEQDB_H_
#define _SEQDB_H_

#include "Time.h"
#include "Gram.h"
#include "GramList.h"
#include "CountFilter.h"
#include "Query.h"

#include <iostream>
#include <fstream>

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <stdint.h>
//add the string.h for strcpy in 446
#include <string.h>
#include <string>
#include <queue>

using namespace std;

#ifdef _WIN32
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

class Candi {
public:
  unsigned id;
  unsigned count;

  Candi(unsigned id, unsigned count) : id(id), count(count) {}
  ~Candi(){}
};

struct queue_entry{
	unsigned m_sid;
	unsigned m_dist;

	queue_entry(unsigned sid, unsigned dist){ m_sid = sid; m_dist = dist; };

	bool operator < (const queue_entry &other) const
	{
		return (m_dist < other.m_dist);
	};
};

struct min_queue_entry{
	unsigned m_sid;
	unsigned m_dist;

	min_queue_entry(unsigned sid, unsigned dist){ m_sid = sid; m_dist = dist; };

	bool operator< (const min_queue_entry &other) const
	{
		return (m_dist>other.m_dist);
	};
};

template <class InvList>
class CSeqDB;

template <class InvList = vector< pair<unsigned, int> > >
class CSeqDB
{
protected:
#ifdef _WIN32
  typedef typename unordered_map<unsigned, CGramList<InvList>* > GramList;
  typedef typename unordered_map<unsigned, string> GramIndex;
  typedef typename unordered_map<unsigned, unsigned> FQueue;
#else
  typedef typename tr1::unordered_map<unsigned, CGramList<InvList>* > GramList;
  typedef typename tr1::unordered_map<unsigned, string> GramIndex;
  typedef typename tr1::unordered_map<unsigned, unsigned> FQueue;
#endif

	vector<string>* data;
	vector<queue_entry> datasizes;
	int* datasize_groupindex;

	GramList gramListUpper;
	vector<unsigned> gramCodeUpper;
	vector<string> gramStrUpper;
	GramList gramListLower;

	int** edcost;
	bool* processedData;
	unsigned** dataCount;
	unsigned fk;

	vector<int> post_cand_low;
	vector<int> post_cand_up;
	// For statistics
public:
	priority_queue<queue_entry> m_queue;
	CQuery* theQuery;
	CGram* gramGenUpper;
	CGram* gramGenLower;
	CCountFilter *filter;

	unsigned gramGenUpper_len;
	unsigned gramGenLower_len;
	bool stop;
	unsigned max_listlen;
	unsigned min_listlen;
	unsigned processed;
	unsigned data_maxlen;
	unsigned gram_maxed;

	//compare with the processed for verification
	unsigned sizeofcandis;

public:
    CSeqDB(const unsigned& maxlen, const unsigned max_gramed, vector<string>* sequences = NULL, CGram* gramgenupper = NULL, CGram* gramgenlower = NULL, CCountFilter *fltable = NULL) :
	  data(sequences), gramGenUpper(gramgenupper), gramGenLower(gramgenlower), filter(fltable)
	{
		this->datasize_groupindex = new int[maxlen+1];
		this->processedData = new bool[sequences->size()];
		for (unsigned i = 0; i < sequences->size(); i++) {
			datasizes.push_back(queue_entry(i, sequences->at(i).length()));
			this->datasize_groupindex[sequences->at(i).length()] = i;
			this->processedData[i] = false;
		}
		this->processed = 0;
		this->gram_maxed = max_gramed;
		this->dataCount = NULL;
		this->data_maxlen = maxlen+1;
		this->edcost = new int *[this->data_maxlen];
		for (unsigned i = 0; i < this->data_maxlen; i++)
			this->edcost[i] = new int[this->data_maxlen];

        //codes for CSeqDB object initial with gramgenupper and gramgenlower parameters
        if (gramGenUpper != NULL){
            this->gramGenUpper_len = this->gramGenUpper->getGramLength();
            this->gramGenLower_len = this->gramGenUpper->getGramLength();
            this->dataCount = new unsigned *[this->gram_maxed];
            for (unsigned i = 0; i < this->gram_maxed; i++) {
                this->dataCount[i] = new unsigned[this->data->size()];
                for (unsigned j = 0; j < this->data->size(); j++)
                    this->dataCount[i][j] = 0;
            }
        }
	}
	~CSeqDB(void)
	{
		for (unsigned i = 0; i < this->data_maxlen; i++)
		{
			if (this->edcost[i] != NULL)
			{
				delete[] this->edcost[i];
				this->edcost[i] = NULL;
			}
		}
		delete[] this->edcost;
		this->edcost = NULL;

		if (this->dataCount != NULL) {
			for (unsigned i = 0; i < this->gram_maxed; i++) {
				if (this->dataCount[i] != NULL)
				{
					delete[] this->dataCount[i];
					this->dataCount[i] = NULL;
				}
			}
			delete[] this->dataCount;
			this->dataCount = NULL;
		}

		if (this->processedData != NULL) {
			delete[] this->processedData;
			this->processedData = NULL;
		}

		if (this->datasize_groupindex != NULL) {
			delete[] this->datasize_groupindex;
			this->datasize_groupindex = NULL;
		}

		this->data = NULL;
		this->gramGenUpper = NULL;
		this->gramGenLower = NULL;
		this->filter = NULL;
	}

	void initParas(const unsigned& kf, CGram* gramgenupper, CGram* gramgenlower, CCountFilter *fltable)
	{
		this->fk = kf;
		this->gramGenUpper = gramgenupper;
		this->gramGenLower = gramgenlower;
		this->filter = fltable;
	}

	void reset()
	{
		for (unsigned i = 0; i < this->data->size(); i++) {
			this->processedData[i] = false;
			for (unsigned j = 0; j < this->gram_maxed; j++)
				this->dataCount[j][i] = 0;
		}
		if (!this->m_queue.empty()) {
			priority_queue<queue_entry> tmp;
			swap(this->m_queue, tmp);
		}
		if (this->post_cand_low.size() > 0) {
			vector<int> tmp;
			swap(this->post_cand_low, tmp);
		}
		if (this->post_cand_up.size() > 0) {
			vector<int> tmp;
			swap(this->post_cand_up, tmp);
		}
		this->processed = 0;
		//this->stop = false;
	}

	void buildindex();
	void insertSequenceIntoIndex(const string& current, const unsigned sid, GramIndex& gramIndex);
	void insertNgramIntoIndex(const string& current, const unsigned gid);

	void printindex();

	bool getEditDistance(const string &tString, const string &qString, int constraint);
	int getRealEditDistance_ns(const string &tString, const string &qString, int constraint);
	int getRealEditDistance_nsd(const string &tString, const string &qString);

	void getGramLists(const vector<unsigned>& gramCodes, GramList& theGramLists, vector<InvList*>& invLists);
	void getGramLists(const vector< pair<unsigned, int> >& gramCodes, GramList& theGramList, vector<InvList*>& invLists);

	static bool cmpInvList(const InvList* a, const InvList* b){
		return a->size() < b->size();
	}
	static bool cmpPosition(const pair<unsigned, int>& a, const pair<unsigned, int>& b) {
		return a.second < b.second;
	}
	void expProbe(const typename InvList::iterator start,
		const typename InvList::iterator end,
		typename InvList::iterator& lbound,
		typename InvList::iterator& ubound,
		unsigned value) const;
	void merge(vector<InvList*>& invLists, const unsigned threshold, vector<Candi>& candiscur);
	void gramQuery(const CQuery& query, const unsigned ged, vector< pair<int, vector<unsigned> > >& similarGrams);
	void getIDBounds(const int querysize, const int threshold, int& idlow, int& idhigh);

	// The fuction for pipe knn search
	void init_threshold(CQuery& query);
	void accumulateFrequency(long ged);
	void knn_postprocess();

	void old_version_knn_postprocess();
};

template <class InvList>
void CSeqDB<InvList>::buildindex()
{
	GramIndex upperGrams;
	for(unsigned sid = 0; sid < this->data->size(); sid++) {
		string current = this->data->at(sid);
		insertSequenceIntoIndex(current, sid, upperGrams);
	}

	//sort the inverted list by position for each upperGram
	typename GramList::iterator listIter;
	for (listIter = this->gramListUpper.begin(); listIter != gramListUpper.end(); listIter++)
		sort(listIter->second->getArray()->begin(), listIter->second->getArray()->end(), cmpPosition);

	typename GramIndex::iterator indexIter;
	unsigned i = 0;

	for(indexIter = upperGrams.begin(); indexIter != upperGrams.end(); indexIter++, i++)
	{
		this->gramCodeUpper.push_back(indexIter->first);
		this->gramStrUpper.push_back(indexIter->second);
		insertNgramIntoIndex(indexIter->second, i);
	}

	//sort the invered list by position for each lowerGram
	for (listIter = this->gramListLower.begin(); listIter != gramListLower.end(); listIter++)
		sort(listIter->second->getArray()->begin(), listIter->second->getArray()->end(), cmpPosition);
}

template <class InvList>
void CSeqDB<InvList>::insertSequenceIntoIndex(const string& current, const unsigned sid, GramIndex& uppergrams)
{
	// add sid to gramListUpper
	vector<string> grams;
	vector< pair<unsigned, int> > gramCodes;
	this->gramGenUpper->decompose(current, grams, gramCodes);
	for(unsigned i = 0; i < gramCodes.size(); i++) {
		pair<unsigned, int> gramCode = gramCodes[i];

		if (this->gramListUpper.find(gramCode.first) == this->gramListUpper.end()) {
			// a new gram
			CGramList<InvList>* newGramList = new CGramList<InvList>();
			this->gramListUpper[gramCode.first] = newGramList;
			newGramList->getArray()->push_back(make_pair(sid, gramCode.second));
		}
		else { 
			// an existing gram
			CGramList<InvList>* theGramList = this->gramListUpper[gramCode.first];
			theGramList->getArray()->push_back(make_pair(sid, gramCode.second));
		}

		if (uppergrams.find(gramCode.first) == uppergrams.end()) {
			uppergrams[gramCode.first] = grams[i];
		}
	}
}

template <class InvList>
void CSeqDB<InvList>::insertNgramIntoIndex(const string& current, const unsigned gid)
{
	// add sid to gramListUpper
	vector< pair<unsigned, int> > gramCodes;
	this->gramGenLower->decompose(current, gramCodes);
	for(unsigned i = 0; i < gramCodes.size(); i++) {
		pair<unsigned, int> gramCode = gramCodes[i];

		if (this->gramListLower.find(gramCode.first) == this->gramListLower.end()) {
			// a new gram
			CGramList<InvList>* newGramList = new CGramList<InvList>();
			this->gramListLower[gramCode.first] = newGramList;
			newGramList->getArray()->push_back(make_pair(gid, gramCode.second));
		}
		else {
			// an existing gram
			CGramList<InvList>* theGramList = this->gramListLower[gramCode.first];
			theGramList->getArray()->push_back(make_pair(gid, gramCode.second));
		}
	}
}

template <class InvList>
void CSeqDB<InvList>::printindex() {
	//print the gramListUpper
	cout << "The upper-level index: " << endl;
	typename GramList::iterator iter;
	for (iter = gramListUpper.begin(); iter != gramListUpper.end(); iter++) {
		cout << iter->first << ":" << endl;
		InvList* theList = iter->second->getArray();
		for (unsigned i = 0; i < theList->size(); i++)
			cout << "(" << theList->at(i).first << "," << theList->at(i).second << ")" << "  ";
		cout << endl;
	}

	//print the gramListLower
	//wait...
}

// Origin verify function from Flamingo
template <class InvList>
bool CSeqDB<InvList>::getEditDistance(const string &tString, const string &qString, int constraint)
{
	int qStrLen = (int)qString.length();
	int tStrLen = (int)tString.length();
	int i, j, min_cost;
	for(i = 0; i < qStrLen+1; i++)
		edcost[0][i] = i;
	int startPos, endPos;
	for(i = 1; i < tStrLen+1; i++) {
		edcost[i][0] = i;
		min_cost = qStrLen;
		startPos = (i - (constraint + 1)) > 1 ? (i - (constraint + 1)) : 1;
		endPos = (i + (constraint + 1)) < qStrLen ? (i + (constraint + 1)) : qStrLen;
		if(startPos > 1) {
			edcost[i][startPos-1] = edcost[i-1][startPos-1] < edcost[i-1][startPos-2] ? (edcost[i-1][startPos-1] + 1) : (edcost[i-1][startPos-2] + 1);
		}
		for(j = startPos; j < endPos+1; j++) {
			edcost[i][j] = std::min(edcost[i-1][j-1]+((tString[i-1]==qString[j-1])?0:1), std::min(edcost[i-1][j]+1, edcost[i][j-1]+1));
			min_cost = edcost[i][j] < min_cost ? edcost[i][j] : min_cost;
		}
		if(j < qStrLen+1) {
			edcost[i][j] = edcost[i][j-1]<edcost[i-1][j-1]?(edcost[i][j-1]+1):(edcost[i-1][j-1]+1);
		}

		if(min_cost > constraint)
			return false;
	}
	return edcost[tStrLen][qStrLen] <= constraint;
}

// Origin verify function from Flamingo
template <class InvList>
int CSeqDB<InvList>::getRealEditDistance_ns(const string &tString, const string &qString, int constraint)
{
	const int qStrLen=(int)qString.length(), tStrLen=(int)tString.length();
	int i, j, min_cost;
	for(i = 0; i < qStrLen+1; i++)
		edcost[0][i] = i;
	int startPos, endPos;
	for(i = 1; i < tStrLen+1; i++) {
		edcost[i][0] = i;
		min_cost = qStrLen;
		startPos = (i - (constraint + 1)) > 1 ? (i - (constraint + 1)) : 1;
		endPos = (i + (constraint + 1)) < qStrLen ? (i + (constraint + 1)) : qStrLen;
		if(startPos > 1) {
			edcost[i][startPos-1] = edcost[i-1][startPos-1] < edcost[i-1][startPos-2] ? (edcost[i-1][startPos-1] + 1) : (edcost[i-1][startPos-2] + 1);
		}
		for(j = startPos; j < endPos+1; j++) {
			edcost[i][j] = std::min(edcost[i-1][j-1]+((tString[i-1]==qString[j-1])?0:1), std::min(edcost[i-1][j]+1, edcost[i][j-1]+1));
			min_cost = edcost[i][j] < min_cost ? edcost[i][j] : min_cost;
		}
		if(j < qStrLen+1) {
			edcost[i][j] = edcost[i][j-1]<edcost[i-1][j-1]?(edcost[i][j-1]+1):(edcost[i-1][j-1]+1);
		}
		if(min_cost > constraint) {
			return min_cost;
		}

	}
	return edcost[tStrLen][qStrLen];
}

//get the real editdistance of two strings
template <class InvList>
int CSeqDB<InvList>::getRealEditDistance_nsd(const string &s1, const string &s2)
{
	const size_t len1 = s1.length(), len2 = s2.length();
	edcost[0][0] = 0;
	for(unsigned int i = 1; i <= len1; ++i) edcost[i][0] = i;
	for(unsigned int i = 1; i <= len2; ++i) edcost[0][i] = i;

	for(unsigned int i = 1; i <= len1; ++i)
		for(unsigned int j = 1; j <= len2; ++j)
                      edcost[i][j] = min(min(edcost[i - 1][j] + 1,edcost[i][j - 1] + 1),
                                          edcost[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) );
	return edcost[len1][len2];
}

template <class InvList>
void CSeqDB<InvList>::getGramLists(const vector<unsigned>& gramCodes, GramList& theGramList, vector<InvList*> &invLists)
{
	unsigned gramCode;
	for (unsigned i = 0; i < gramCodes.size(); i++) {
		gramCode = gramCodes[i];
		if (theGramList.find(gramCode) != theGramList.end()) {
			InvList* tmp = theGramList[gramCode]->getArray();
			invLists.push_back(tmp);
		}
	}
}

template <class InvList>
void CSeqDB<InvList>::getGramLists(const vector< pair<unsigned, int> >& gramCodes, GramList& theGramList, vector<InvList*> &invLists)
{
	unsigned gramCode;
	for(unsigned i = 0; i < gramCodes.size(); i++) {
		gramCode = gramCodes[i].first;
		if(theGramList.find(gramCode) != theGramList.end()) {
		  InvList* tmp = theGramList[gramCode]->getArray();
		  invLists.push_back(tmp);
		}
	}
}

template <class InvList>
void CSeqDB<InvList>::expProbe(const typename InvList::iterator start,
	 const typename InvList::iterator end,
	 typename InvList::iterator& lbound,
	 typename InvList::iterator& ubound,
	 unsigned value) const {

	unsigned c = 0;
	lbound = start;
	for(;;) {
		ubound = start + (1 << c);
		if(ubound >= end) {
			ubound = end;
			return;
		}
		else if((*ubound).first >= value) return;
		c++;
		lbound = ubound;
	}
}

//str which in the invLists and count > threshold, put its code into the candiscur
template <class InvList>
void CSeqDB<InvList>::merge(vector<InvList*>& invLists, const unsigned threshold, vector<Candi>& candiscur)
{
    //size: smaller ~ bigger
	sort(invLists.begin(), invLists.end(), CSeqDB<InvList>::cmpInvList);

	unsigned numShortLists = invLists.size() - threshold + 1;

	// process the short lists using the algorithm from
	// Naoaki Okazaki, Jun'ichi Tsujii
	// "Simple and Efficient Algorithm for Approximate Dictionary Matching", COLING 2010

    //candis for store all the Candi involved in invLists
	vector<Candi> candis;
	for(unsigned i = 0; i < numShortLists; i++) {
		vector<Candi> tmp;

		tmp.reserve(candis.size() + invLists[i]->size());
		vector<Candi>::const_iterator candiIter = candis.begin();

		typename InvList::iterator invListIter = invLists[i]->begin();
		typename InvList::iterator invListEnd = invLists[i]->end();

        //merge candis and invListIter, accumulate the counts
		while(candiIter != candis.end() || invListIter != invListEnd) {
			if(candiIter == candis.end() || (invListIter != invListEnd && candiIter->id > (*invListIter).first)) {
				tmp.push_back(Candi((*invListIter).first, 1));
				invListIter++;
			}
			else if (invListIter == invListEnd || (candiIter != candis.end() && candiIter->id < (*invListIter).first)) {
				tmp.push_back(Candi(candiIter->id, candiIter->count));
				candiIter++;
			}
			else {
				tmp.push_back(Candi(candiIter->id, candiIter->count + 1));
				candiIter++;
				invListIter++;
			}
		}
		std::swap(candis, tmp);
	}

	// process long lists
	unsigned stop = invLists.size();
	for(unsigned j = numShortLists; j < stop; j++) {
		if(candis.size() < stop - j) break; // termination heuristic

		vector<Candi> tmp;
		tmp.reserve(candis.size());
		typename InvList::iterator listIter = invLists[j]->begin();

		for(unsigned i = 0; i < candis.size(); i++) {
			typename InvList::iterator start, end;
			expProbe(listIter, invLists[j]->end(), start, end, candis[i].id);
			//lower_bound
			//Returns an iterator pointing to the first element in the range [first,last)
			//which does not compare less than val
			//listIter = lower_bound(start, end, candis[i].id);
			for (listIter = start; listIter != end; listIter++)
				if ((*listIter).first >= candis[i].id)
					break;
			if(listIter != invLists[j]->end() && (*listIter).first == candis[i].id) {
				candis[i].count++;
				if(candis[i].count >= threshold) candiscur.push_back(candis[i]);
				else tmp.push_back(candis[i]);
			}
			else {
				if(candis[i].count + stop - j > threshold)
					tmp.push_back(candis[i]);
			}
		}

		std::swap(candis, tmp);
	}

	for (unsigned i = 0; i < candis.size(); i++) {
		candiscur.push_back(candis[i]);
	}
}

//query => ngrams => lower str list => merge count => upper str list
template <class InvList>
void CSeqDB<InvList>::gramQuery(const CQuery& query, const unsigned ged, vector< pair<int, vector<unsigned> > >& similarGrams)
{
	const vector<string>& gramStrs = query.ngrams;
	vector< pair<unsigned, int> > gramCodes;
	//class Candi: contains its id and its counts
	vector<Candi> candis;
	vector<unsigned> results;
	int threshold = (int)this->gramGenLower->getGramLength();
	threshold = (int)(this->gramGenUpper->getGramLength()) - threshold + 1 - (int)ged*threshold;
	for (unsigned i = 0; i < gramStrs.size(); i++) {
		gramCodes.clear();
		this->gramGenLower->decompose(gramStrs[i], gramCodes);
		vector<InvList*> lists;
		this->getGramLists(gramCodes, this->gramListLower, lists);
		results.clear();
		if ((unsigned)threshold > lists.size()) {
			continue;
		}

		candis.clear();
		//merge the lists
		this->merge(lists, (unsigned)threshold, candis);
		for(unsigned j = 0; j < candis.size(); j++) {
			if (this->getEditDistance(this->gramStrUpper[candis[j].id], gramStrs[i], (int)ged))
				results.push_back(this->gramCodeUpper[candis[j].id]);
		}
		similarGrams.push_back(make_pair(this->theQuery->gramCodes[i].second, results));
	}
}

template <class InvList>
void CSeqDB<InvList>::init_threshold(CQuery& query) {
	int idlow, idhigh, ed;
	this->processed = 0;
	this->getIDBounds((int)query.length, 0, idlow, idhigh);
	idhigh = idhigh >(idlow + (int)query.constraint) ? idhigh : ((idlow + (int)query.constraint) < (int)this->data->size() ? (idlow + (int)query.constraint) : (int)this->data->size());
	idlow = idlow < (idhigh - (int)query.constraint) ? idlow : ((idhigh - (int)query.constraint)>0 ? (idhigh - (int)query.constraint) : 0);
	for (; idlow <= idhigh; idlow++) {
		if (this->m_queue.size() < query.constraint) {
			ed = this->getRealEditDistance_nsd(this->data->at(idlow), query.sequence);
			this->processedData[idlow] = true;
			this->processed++;
			this->m_queue.push(queue_entry(idlow, (unsigned)ed));
			query.threshold = this->m_queue.top().m_dist;
			continue;
		}
		//new modify for simple index, other index can use it too.
		//where there are enough (query.constraint) candidates in the m_queue for the threshold computing
		//if (this->m_queue.size() == query.constraint)
        //    return ;

		if (query.threshold == 0) {
			return;
		}

		ed = this->getRealEditDistance_ns(this->data->at(idlow), query.sequence, (int)query.threshold);
		this->processedData[idlow] = true;
		this->processed++;
		if (ed < (int)query.threshold) {
			this->m_queue.pop();
			this->m_queue.push(queue_entry(idlow, (unsigned)ed));
			query.threshold = this->m_queue.top().m_dist;
			if (query.threshold == 0) {
				return;
			}
		}
	}
	query.threshold = this->m_queue.top().m_dist;
}

//knn_postprocess, older version, until 2015-11-25
template <class InvList>
void CSeqDB<InvList>::old_version_knn_postprocess() {
	int ed, idlow, idhigh;
	//unsigned j;
	idlow = idhigh = -1;
	this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idlow, idhigh);

	if (idlow <= idhigh) {
		for (int dataCode = idlow; dataCode <= idhigh; dataCode++) {
			if (this->processedData[dataCode])
				continue;

			if (this->m_queue.size() < this->theQuery->constraint) {
				ed = this->getRealEditDistance_nsd(this->data->at(dataCode), this->theQuery->sequence);
				this->processedData[dataCode] = true;
				this->processed++;
				this->m_queue.push(queue_entry(dataCode, (unsigned)ed));
				this->theQuery->threshold = this->m_queue.top().m_dist;
				continue;
			}

			// Check the bounds on approximate ngrams
			if (!this->processedData[dataCode] && (int)this->dataCount[0][dataCode] >= this->filter->tabUpQuery[this->theQuery->threshold][0]) {
				ed = this->getRealEditDistance_ns(this->data->at(dataCode), this->theQuery->sequence, (int)this->theQuery->threshold);
				this->processedData[dataCode] = true;
				this->processed++;
				if (ed < (int)this->theQuery->threshold) {
					if (this->m_queue.size() == this->theQuery->constraint)
						this->m_queue.pop();
					queue_entry tmp_entry(dataCode, (unsigned)ed);
					this->m_queue.push(tmp_entry);
					this->theQuery->threshold = this->m_queue.top().m_dist;
					if (this->theQuery->threshold == 0) {
						//this->stop = true;
						return;
					}

					this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, dataCode, idhigh);
				}
			}
		}
	}
}

template <class InvList>
void CSeqDB<InvList>::knn_postprocess() {
	bool checked;
	unsigned i, j;
	int ed, appgram_count, idlow, idhigh, idnew;
	this->getIDBounds((int)this->theQuery->length, this->theQuery->threshold, idlow, idhigh);
	const int& querylen_index = this->datasize_groupindex[this->theQuery->length];
	if (idlow <= idhigh) {
		// For those sizes no larger than query
		for (; idlow <= querylen_index; idlow++) {
			if (!this->processedData[idlow]) {
				if (this->dataCount[0][idlow] == 0) {
					// Even do not have a similar gram
					if (this->dataCount[1][idlow] == 0)
						continue;
					//approximate ngrams
					this->post_cand_low.push_back(idlow);
					continue;
				}
				// Check the bounds on approximate ngrams
				checked = true;
				appgram_count = 0;
				for (i = 0; i < this->gram_maxed; i++) {
					appgram_count += this->dataCount[i][idlow];
					if (appgram_count < this->filter->tabUpQuery[this->theQuery->threshold][i]) {
						checked = false;
						break;
					}
				}
				if (checked) {
					ed = this->getRealEditDistance_ns(this->data->at(idlow), this->theQuery->sequence, (int)this->theQuery->threshold);
					this->processedData[idlow] = true;
					this->processed++;
					if (ed < (int)this->theQuery->threshold) {
						this->m_queue.pop();
						this->m_queue.push(queue_entry(idlow, (unsigned)ed));
						if (this->m_queue.top().m_dist < this->theQuery->threshold) {
							this->theQuery->threshold = this->m_queue.top().m_dist;
							if (this->theQuery->threshold == 0)
								return;
							this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idnew, idhigh);
							if (idnew > querylen_index) {
								idlow = idnew;
								break;
							}
							idlow = idlow < (idnew - 1) ? (idnew - 1) : idlow;
						}
					}
				}
			}
		}
		// For those sizes larger than query
		for (; idlow <= idhigh; idlow++) {
			if (!this->processedData[idlow]) {
				if (this->dataCount[0][idlow] == 0) {
					// Even do not have a similar gram
					if (this->dataCount[1][idlow] == 0)
						continue;
					this->post_cand_up.push_back(idlow);
					continue;
				}
				// Check the bounds on approximate ngrams
				checked = true;
				appgram_count = 0;
				const int& tmp_datasize = this->data->at(idlow).size();
				for (i = 0; i < this->gram_maxed; i++) {
					appgram_count += this->dataCount[i][idlow];
					if ((appgram_count + this->filter->tabUp[this->theQuery->threshold][i]) < tmp_datasize) {
						checked = false;
						break;
					}
				}
				if (checked) {
					ed = this->getRealEditDistance_ns(this->data->at(idlow), this->theQuery->sequence, (int)this->theQuery->threshold);
					this->processedData[idlow] = true;
					this->processed++;
					if (ed < (int)this->theQuery->threshold) {
						this->m_queue.pop();
						this->m_queue.push(queue_entry(idlow, (unsigned)ed));
						if (this->m_queue.top().m_dist < this->theQuery->threshold) {
							this->theQuery->threshold = this->m_queue.top().m_dist;
							if (this->theQuery->threshold == 0)
								return;
							this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idnew, idhigh);
							if (idnew >= idhigh)
								break;
							idlow = idlow < (idnew - 1) ? (idnew - 1) : idlow;
						}
					}
				}
			}
		}
		// Check approximate ngrams
		this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idlow, idhigh);
		for (j = 0; j < this->post_cand_low.size(); j++) {
			const int& dataid = this->post_cand_low[j];
			if (dataid >= idlow && dataid <= idhigh) {
				checked = true;
				appgram_count = 0;
				for (i = 1; i < this->gram_maxed; i++) {
					appgram_count += this->dataCount[i][dataid];
					if (appgram_count < this->filter->tabUpQuery[this->theQuery->threshold][i]) {
						checked = false;
						break;
					}
				}
				if (checked) {
					ed = this->getRealEditDistance_ns(this->data->at(dataid), this->theQuery->sequence, (int)this->theQuery->threshold);
					this->processedData[dataid] = true;
					this->processed++;
					if (ed < (int)this->theQuery->threshold) {
						this->m_queue.pop();
						this->m_queue.push(queue_entry(dataid, (unsigned)ed));
						if (this->theQuery->threshold > this->m_queue.top().m_dist) {
							this->theQuery->threshold = this->m_queue.top().m_dist;
							if (this->theQuery->threshold == 0) {
								return;
							}
							this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idlow, idhigh);
						}
					}
				}
			}
		}
		for (j = 0; j < this->post_cand_up.size(); j++) {
			const int& dataid = this->post_cand_up[i];
			if (dataid >= idlow && dataid <= idhigh) {
				checked = true;
				appgram_count = 0;
				const int& tmp_datasize = (int)this->data->at(dataid).size();
				for (i = 1; i < this->gram_maxed; i++) {
					appgram_count += this->dataCount[i][dataid];
					if ((appgram_count + this->filter->tabUp[this->theQuery->threshold][i]) < tmp_datasize) {
						checked = false;
						break;
					}
				}
				if (checked) {
					ed = this->getRealEditDistance_ns(this->data->at(dataid), this->theQuery->sequence, (int)this->theQuery->threshold);
					this->processedData[dataid] = true;
					this->processed++;
					if (ed < (int)this->theQuery->threshold) {
						this->m_queue.pop();
						this->m_queue.push(queue_entry(dataid, (unsigned)ed));
						if (this->theQuery->threshold > this->m_queue.top().m_dist) {
							this->theQuery->threshold = this->m_queue.top().m_dist;
							if (this->theQuery->threshold == 0) {
								return;
							}
							this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idlow, idhigh);
						}
					}
				}
			}
		}
	}
}

// Accumulate the frequency of approximate ngrams
// the largest value of ged is decided by upper gram length according to the CA strategy
template <class InvList>
void CSeqDB<InvList>::accumulateFrequency(long ged)
{
	typename InvList::iterator iterlow, iterhigh;

	if (ged == 0) {
		vector<InvList*> lists;
		this->getGramLists(this->theQuery->gramCodes, this->gramListUpper, lists);

		for (unsigned i = 0; i < lists.size(); i++) {
			iterlow = lists[i]->begin();
			iterhigh = lists[i]->end();
			for (; iterlow != iterhigh; iterlow++)
				if (abs(this->theQuery->gramCodes[i].second - (*iterlow).second) <= this->theQuery->threshold)
					this->dataCount[ged][(*iterlow).first]++;
		}
	}
	else {
		vector<InvList*> lists;
		vector< pair<int, vector<unsigned> > > similarGrams;
		this->gramQuery(*(this->theQuery), ged, similarGrams);
		
		for (unsigned i = 0; i < similarGrams.size(); i++) {
			vector<unsigned> v(similarGrams[i].second);
			this->getGramLists(v, this->gramListUpper, lists);

			for (unsigned j = 0; j < lists.size(); j++) {
				iterlow = lists[j]->begin();
				iterhigh = lists[j]->end();
				for (; iterlow != iterhigh; iterlow++)
					if (abs(similarGrams[i].first - (*iterlow).second) <= this->theQuery->threshold)
						this->dataCount[ged][(*iterlow).first]++;
			}
		}
	}
}

template <class InvList>
void CSeqDB<InvList>::getIDBounds(const int querysize, const int threshold, int& idlow, int& idhigh)
{
	int lowlen = querysize - threshold > 0 ? querysize - threshold : 0;
	int highlen = querysize + threshold;
	if ((unsigned)lowlen > this->datasizes.front().m_dist) {
		vector<queue_entry>::iterator _iterLow = lower_bound(this->datasizes.begin(), this->datasizes.end(), queue_entry(0, (unsigned)lowlen));
		idlow = (int)((*_iterLow).m_sid);
	}
	else {
		idlow = this->datasizes.front().m_sid;
	}
	if ((unsigned)highlen < this->datasizes.back().m_dist) {
		vector<queue_entry>::iterator _iterHigh = upper_bound(this->datasizes.begin(), this->datasizes.end(), queue_entry(0, (unsigned)highlen));
		idhigh = (int)((*_iterHigh).m_sid);
	}else {
		idhigh = this->datasizes.back().m_sid;
	}
}

#endif /* _SEQDB_H_ */
