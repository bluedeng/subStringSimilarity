#include "Time.h"
#include "Gram.h"
#include "CountFilter.h"
#include "SeqDB.h"
#include "Query.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>
#include <set>
#include <string>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

using namespace std;

int gl_qmax = 5;
int gl_qmin = 2;
int gl_topk = 1;
int gl_topk_f = 1;
int gl_threshold;
bool output = false;
vector<string> dataset;
vector<string> queryset;
unsigned gl_maxLen = 0;

//double queryPrep_time = 0.0;
double index_time = 0.0;
double acc_time = 0.0;
double query_time = 0.0;
double all_time = 0.0;
long processednum = 0;

void usage()
{
	fprintf(stdout, "**********************************************\n");
	fprintf(stdout, "SSEARCH 1.0 Copyright 2010, NUS, Singapore\n");
	fprintf(stdout, "Sequence search algorithms based on CA\n");
	fprintf(stdout, "----------------------------------------------\n\n");
	fprintf(stdout, "Usage: topkSubStringMatching [OPTION] ${INPUTDB} ${QUERY} [ed]\n");
	fprintf(stdout, "-o\n   ouput the result set of sequence ids.\n");
	fprintf(stdout, "We will process the topK for 5, 10, 20, 30, 40 \n");
	fprintf(stdout, "\n----------------------------------------------\n");
	fprintf(stdout, "Author : Xiaoli Wang, SOC, NUS, Singapore.\n");
	fprintf(stdout, "GMail  : xiaolistue@gmail.com");
	fprintf(stdout, "\n**********************************************\n");
}

//compare the strings between psearch and dataset, and get some possible arguments
void parseOptions(int argc, const char* argv[])
{
	if (argc > 3) {
		for (int i = 1; argc - 2>i; i++) {
			//compare first size chars of string1 and string2, if true return 0
			if (strncmp(argv[i], "-k", 2) == 0) {
				if (1 != sscanf(argv[++i], "%d", &gl_topk)) {
					fprintf(stdout, "Error occured while parsing the value of k.\n");
					exit(-1);
				}
				if (gl_topk <= 0) {
					fprintf(stdout, "The top k value error : %d.\n", gl_topk);
					exit(-1);
				}
			}

			if (strncmp(argv[i], "-f", 2) == 0) {
				if (1 != sscanf(argv[++i], "%d", &gl_topk_f)) {
					fprintf(stdout, "Error occured while parsing the value of k for f-queue.\n");
					exit(-1);
				}
				if (gl_topk_f <= 0) {
					fprintf(stdout, "The queue k value error : %d.\n", gl_topk_f);
					exit(-1);
				}
			}

			if (strncmp(argv[i], "-o", 2) == 0) {
				output = true;
			}
		}
	}
}

//string : file(strings) -> data
void read(const char *fileName, vector<string> &data)
{
	ifstream _fs;
	_fs.open(fileName, ios::in);
	//failed to read the file
	if (_fs.fail()) {
		fprintf(stderr, "Error: Failed to open the data file - %s\n", fileName);
		fprintf(stderr, "Please check the file name and restart this program.\n");
		fprintf(stderr, "The program will be terminated!\n");
		exit(-1);
	}

	string line;
	while (!_fs.eof()) {
		getline(_fs, line);

		if (line.empty())
			continue;

		for (unsigned i = 0; line.size() > i; i++)
			if (line[i] >= 'A' && line[i] <= 'Z')
				line[i] += 32;

		if (line[line.size() - 1] == '\r')
			line.erase(line.size() - 1);

		if (gl_maxLen < line.size())
			gl_maxLen = line.size();

		data.push_back(line);
	}

	_fs.close();
}

static CSeqDB<> *pThis;
#define PTHEAD_NUM 2

volatile int threadnum = 0;
pthread_mutex_t running_mutex = PTHREAD_MUTEX_INITIALIZER;

// The fuction of the thread
void *thr_fn_knn(void *arg)
{
	//pthread_self() : use it for the thread id
	int iret = pthread_detach(pthread_self());
	if (iret != 0)
	{
		fprintf(stderr, "Can't detach at my thread!\n");
		pthread_mutex_lock(&running_mutex);
		threadnum--;
		pthread_mutex_unlock(&running_mutex);
		pthread_exit(NULL);   // Exit the thread
		return ((void *)0);
	}

	// Accumulate frequency
	long ged = (long)arg;
	pThis->accmulateFrequencyForSubString(ged);

	pthread_mutex_lock(&running_mutex);
	threadnum--;
	pthread_mutex_unlock(&running_mutex);
	pthread_exit(NULL);
	return ((void *)0);
}

int pipe_knn(CQuery* query) {
	pThis->theQuery = query;

	class TSINGHUA_CLIPSE_UTIL::TimeRecorder _time;

	pthread_t threads[PTHEAD_NUM];
	int error = 0;
	for (int i = 0; i < PTHEAD_NUM; i++)
	{

		pthread_mutex_lock(&running_mutex);
		threadnum++;
		pthread_mutex_unlock(&running_mutex);
		error = pthread_create(&threads[i], NULL, thr_fn_knn, (void *)i);
		if (0 != error)
		{
			fprintf(stderr, "Cannot create thread:  %d!\n", error);
			return 0;
		}

	}

	// Check if all threads are closed.
	while (true)
	{
		if (threadnum <= 0)
		{
			_time.check();
			acc_time += _time.diffTime(0, 1);

			// Try post precessing here with more approximate bounds
			class TSINGHUA_CLIPSE_UTIL::TimeRecorder _time1;
			pThis->subString_range(gl_threshold);
			_time1.check();
			query_time += _time1.diffTime(0, 1);
			// Then end of whole processing here
			return 1;
		}
		sleep(0.0001);
	}
	return 0;
}

//change int to string
string changeI2S(int i) {
	stringstream ss;
	string temp;
	ss << i;
	ss >> temp;
	ss.clear();
	return temp;
}

//change string to int
int changeS2I(string temp) {
	stringstream ss;
	int i;
	ss << temp;
	ss >> i;
	ss.clear();
	return i;
}

int main(int argc, const char* argv[]) {
	//./rangeSubStringMatch dataset query ed
	if (argc < 4) {
		usage();
		exit(-1);
	}

	parseOptions(argc, argv);

	// read data file
	read(argv[argc - 3], dataset);
	if (dataset.size() == 0) {
		fprintf(stderr, "Error: Failed to open the dataset file - %s\n", argv[argc - 2]);
		fprintf(stderr, "Please check the file name and restart this program.\n");
		fprintf(stderr, "The program will be terminated!\n");
		exit(-1);
	}

	// read query file
	read(argv[argc - 2], queryset);
	if (queryset.size() == 0) {
		fprintf(stderr, "Error: Failed to open the query file - %s\n", argv[argc - 1]);
		fprintf(stderr, "Please check the file name and restart this program.\n");
		fprintf(stderr, "The program will be terminated!\n");
		exit(-1);
	}

	gl_threshold = changeS2I((string)argv[argc - 1]);
	string stemp = argv[argc - 2];
	ofstream fout((stemp + "_subString_range").c_str());
	fout << "threshold, Average_candidate_number, Average_indextime, Average_acctime, Average_querytime, Average_time" << endl;

	//in case of we build the index and use it in the same program
	CGram gramgenupper(gl_qmax, false);
	CGram gramgenlower(gl_qmin, false);
	CCountFilter thefilter(&gramgenupper, gl_maxLen, PTHEAD_NUM);
	CSeqDB<> db_load(gl_maxLen, PTHEAD_NUM, &dataset, &gramgenupper, &gramgenlower, &thefilter);
	class TSINGHUA_CLIPSE_UTIL::TimeRecorder _time1;
	//build the index and time it
	db_load.buildindex();
	_time1.check();
	index_time += _time1.diffTime(0, 1);

	pThis = &db_load;

	ofstream fdetail((stemp + "_subString_range_" + changeI2S(gl_threshold)).c_str());

	acc_time = 0.0;
	query_time = 0.0;
	all_time = 0.0;
	processednum = 0;

	for (unsigned i = 0; i < queryset.size(); i++) {//get each query
		if (queryset[i].size() < db_load.gramGenUpper_len) {
			fprintf(stdout, "We support queries with length larger than %d\n", db_load.gramGenUpper_len);
			fprintf(stdout, "Too small query size for [%s]\n", queryset[i].c_str());
			continue;
		}
		//cout the processing query index
		cout << "processing the " << i << "th" << endl;
		db_load.reset();//reset the db_load

		CQuery query(queryset[i].c_str(), db_load.gramGenUpper, (unsigned)gl_topk);
		class TSINGHUA_CLIPSE_UTIL::TimeRecorder _time;
		query.threshold = gl_threshold;
		db_load.filter->setQueryLb(query);
		pipe_knn(&query);
		_time.check();
		all_time += _time.diffTime(0, 1);
		processednum += db_load.processed;

		fdetail << "Query for: " << queryset[i] << endl;
		if (output) {
			fprintf(stdout, "Finish the query for [%s]:\n", queryset[i].c_str());
			while (!db_load.m_queue.empty()) {
				const queue_entry &entry = db_load.m_queue.top();
				fprintf(stdout, "%d: [%s]\t%d\n", entry.m_sid, dataset[entry.m_sid].c_str(), entry.m_dist);
				fdetail << dataset[entry.m_sid].c_str() << " " << entry.m_dist << endl;
				db_load.m_queue.pop();
			}
			fprintf(stdout, "**********************************************\n");
		}
	}
	fdetail.close();

	//output
	fout << gl_threshold  << "  " << processednum / queryset.size() << "  " << index_time << "  " 
		 << acc_time / queryset.size() << "  " << query_time / queryset.size() << "  " << all_time / queryset.size() << endl;

	fout.close();
	return 0;
}