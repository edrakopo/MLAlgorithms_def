#include <algorithm>
#include <vector>
#include <numeric>
#include <utility>
#include <iostream>
#include <iterator>
#include <list>
#include "TObject.h"

/* 
 std::vector<double> values;
 MakeBins thebins = MakeBins(values,5);
 std::cout << thebins.binStart(i) << " " << thebins.binStop(i) << std::endl;
*/

class MakeBins {

 public:

  //  ClassDef(MakeBins,1);

  MakeBins(std::vector<double>& values, const unsigned int nBin);

  MakeBins(std::vector<std::pair< double, double> >& values, std::vector<double>& bins ); 

  MakeBins(std::vector<double>& breaks);


  MakeBins(double beg, double end, int nBin);

  ~MakeBins(){}

  unsigned int toBin(const double value) const;
 
  double binCenter(const unsigned int i) const;

  const std::vector<double> binCenters() const { return m_binCenter; }

  double binStart(const unsigned int i) const;

  double binStop(const unsigned int i) const;

  int size() const;

 private: 

  typedef std::pair<double, double> Entry;
  typedef std::vector<Entry> Entries;
  Entries m_entries;
  std::vector<double> m_binCenter;

};

int MakeBins::size() const {
  return m_entries.size();
}

double MakeBins::binCenter(const unsigned int i) const {
  return m_binCenter[i];
}
 
double MakeBins::binStop(const unsigned int i) const{
  return m_entries[i].second;
}

double MakeBins::binStart(const unsigned int i) const{
  return m_entries[i].first;
}

unsigned int MakeBins::toBin(const double value) const{

  // std::cout << value << std::endl;
  unsigned int entry = 0u;
  for (; entry < m_entries.size(); ++entry){
    if (value > m_entries[entry].first && value < m_entries[entry].second) break;
  }  
  return entry < m_entries.size() ? entry : m_entries.size() -1 ;
}

MakeBins::MakeBins(std::vector<double>& values, const unsigned int nBin){


  std::sort(values.begin(), values.end() );
  unsigned int entrySize = values.size()/nBin;
  std::vector<double>::iterator start = values.begin();
  std::vector<double>::iterator stop;
  for (unsigned int iEntry = 0; iEntry < nBin; ++iEntry){
    std::cout << iEntry << std::endl;
    stop = start + entrySize;
    double mean = std::accumulate(start,stop,0.0)/double(entrySize);
    m_binCenter.push_back(mean);
    m_entries.push_back(std::make_pair(*start, *(stop - 1)));
    start = stop;
  } //iEntry
  

    /*
  unsigned int nEntry = svalues.size()/nBin;
  unsigned int startIndex = 0;
  std::cout << nEntry << std::endl;
  for (i = 0; i < nBin; ++i){
    unsigned int stopIndex = startIndex + nEntry; 
    double mean = 0;
    for (int in = startIndex; in < stopIndex; ++in) mean += svalues[in]; 
    mean /= float(nEntry);
    m_binCenter.push_back(mean);
    std::cout << startIndex <<  " " << stopIndex << std::endl;
    m_entries.push_back(std::make_pair(svalues[startIndex], svalues[stopIndex]));
    startIndex = stopIndex;
  } 
    */

  for (unsigned i = 0; i < m_entries.size(); ++i){
    std::cout << m_entries[i].first << " " << m_entries[i].second <<  " "  <<  m_binCenter[i] << std::endl;
  }
}



MakeBins::MakeBins(std::vector<std::pair<double, double> >& values, 
                          std::vector<double>& bins): m_entries(values), 
                                                      m_binCenter(bins) {}
  
MakeBins::MakeBins(std::vector<double>& breaks){

  for (unsigned int i = 1; i < breaks.size(); ++i){
      m_entries.push_back(std::make_pair(breaks[i-1], breaks[i]));
      m_binCenter.push_back(0.5 * (breaks[i-1] + breaks[i]));
  }   
 
  for (unsigned i = 0; i < m_entries.size(); ++i){
    std::cout << m_entries[i].first << " " << m_entries[i].second <<  " "  <<  m_binCenter[i] << std::endl;
  }

}

MakeBins::MakeBins(double beg, double end, int nBin){

  double binSize = (end - beg) / nBin;
  for (int i = 0; i< nBin; ++i){
    double start = beg + i* binSize;
    double stop = start + binSize;
    m_entries.push_back(std::make_pair(start, stop));     
    m_binCenter.push_back(0.5 * (start + stop));
    cout << start << " " << stop << " " << i << std::endl ;  
  }
 
}

