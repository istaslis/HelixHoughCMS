#include "MLoVetere/HTTrackSeeding/interface/HelixHough.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <sys/time.h>

using namespace std;
using namespace PHENIXHough;


static inline void in_place_counting_sort_unique(vector<unsigned int>& A, vector<unsigned int>& C, unsigned int MAX)
{
  unsigned int SIZE = A.size();
  for(unsigned int i=0;i<SIZE;++i)
  {
    ++C[A[i]];
  }
  unsigned int current = 0;
  for (unsigned int i=0;i<MAX;++i)
  {
    A[current] = ( (((C[i]!=0)-1)&(A[current])) ^ (((C[i]==0)-1)&i) );
    current += (C[i]!=0);
  }
  A.resize(current);
}


static inline void in_place_counting_unique(vector<unsigned int>& A, vector<unsigned int>& C)
{
  unsigned int SIZE = A.size();
  unsigned int current = 0;
  for(unsigned int i=0;i<SIZE;++i)
  {
    C[A[i]] += 1;
    if((C[A[i]]==1))
    {
      A[current] = A[i];
      current += 1;
    }
  }
  A.resize(current);
}


void HelixHough::setRange(const BinEntryPair5D& bp, HelixRange& range1, HelixRange& range2, unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0)
{
  float dzdl_size = (range1.max_dzdl - range1.min_dzdl)/((float)n_dzdl);
  float z0_size = (range1.max_z0 - range1.min_z0)/((float)n_z0);
  float phi_size = (range1.max_phi - range1.min_phi)/((float)n_phi);
  float d_size = (range1.max_d - range1.min_d)/((float)n_d);
  float k_size = (range1.max_k - range1.min_k)/((float)n_k);
  
  unsigned int z0_bin = 0;
  unsigned int dzdl_bin = 0;
  unsigned int k_bin = 0;
  unsigned int d_bin = 0;
  unsigned int phi_bin = 0;
  
  bp.bin5D(n_d, n_k, n_dzdl, n_z0, phi_bin, d_bin, k_bin, dzdl_bin, z0_bin);
  range2.min_phi = range1.min_phi + phi_size*((float)(phi_bin));
  range2.max_phi = range2.min_phi + phi_size;
  range2.min_d = range1.min_d + d_size*((float)(d_bin));
  range2.max_d = range2.min_d + d_size;
  range2.min_k = range1.min_k + k_size*((float)(k_bin));
  range2.max_k = range2.min_k + k_size;
  range2.min_dzdl = range1.min_dzdl + dzdl_size*((float)(dzdl_bin));
  range2.max_dzdl = range2.min_dzdl + dzdl_size;
  range2.min_z0 = range1.min_z0 + z0_size*((float)(z0_bin));
  range2.max_z0 = range2.min_z0 + z0_size;
}


static inline void setClusterRange(HelixRange& range1, HelixRange& range2, ParRange& prange, unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0)
{
  float dzdl_size = (range1.max_dzdl - range1.min_dzdl)/((float)n_dzdl);
  float z0_size = (range1.max_z0 - range1.min_z0)/((float)n_z0);
  float phi_size = (range1.max_phi - range1.min_phi)/((float)n_phi);
  float d_size = (range1.max_d - range1.min_d)/((float)n_d);
  float k_size = (range1.max_k - range1.min_k)/((float)n_k);
  
  range2.min_phi = range1.min_phi + phi_size*((float)(prange.min_phi));
  range2.max_phi = range1.min_phi + phi_size*((float)(prange.max_phi + 1));
  range2.min_d = range1.min_d + d_size*((float)(prange.min_d));
  range2.max_d = range1.min_d + d_size*((float)(prange.max_d + 1));
  range2.min_k = range1.min_k + k_size*((float)(prange.min_k));
  range2.max_k = range1.min_k + k_size*((float)(prange.max_k + 1));
  range2.min_dzdl = range1.min_dzdl + dzdl_size*((float)(prange.min_dzdl));
  range2.max_dzdl = range1.min_dzdl + dzdl_size*((float)(prange.max_dzdl + 1));
  range2.min_z0 = range1.min_z0 + z0_size*((float)(prange.min_z0));
  range2.max_z0 = range1.min_z0 + z0_size*((float)(prange.max_z0 + 1));
}


void HelixHough::findHelices(vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits, vector<SimpleTrack3D>& tracks, unsigned int maxtracks)
{
  vote_time = 0.;
  xy_vote_time = 0.;
  z_vote_time = 0.;
  cluster_time = 0.;
  
  index_mapping.clear();
  index_mapping.resize(hits.size(), 0);
  hit_used->clear();
  for(unsigned int i=0;i<hits.size();i++)
  {
    index_mapping[i] = hits[i].index;
    hits[i].index = i;
  }
  (*hit_used).assign(hits.size(), false);
  
  initEvent(hits, min_hits);
  
  max_tracks = maxtracks;
  base_hits = &hits;
  (*(hits_vec[start_zoom])) = hits;
  current_range = top_range;
  zoomranges.clear();
  for(unsigned int z=0;z<=max_zoom;z++)
  {
    zoomranges.push_back(top_range);
  }
  vector<SimpleTrack3D> temp_tracks;
  if((separate_by_helicity==true) && (only_one_helicity==false))
  {
    helicity=true;
    findHelices(min_hits, max_hits, temp_tracks, maxtracks, start_zoom);
    helicity=false;
    findHelices(min_hits, max_hits, temp_tracks, maxtracks, start_zoom);
  }
  else
  {
    findHelices(min_hits, max_hits, temp_tracks, maxtracks, start_zoom);
  }
  
  for(unsigned int i=0;i<hits.size();i++)
  {
    hits[i].index = index_mapping[i];
  }
  for(unsigned int t=0;t<temp_tracks.size();t++)
  {
    for(unsigned int h=0;h<temp_tracks[t].hits.size();h++)
    {
      if (temp_tracks[t].hits[h].index!=(unsigned)-1)
        temp_tracks[t].hits[h].index = index_mapping[temp_tracks[t].hits[h].index];
    }
  }
  
  finalize(temp_tracks, tracks);
  
  if(print_timings == true)
  {
    cout<<"vote time = "<<vote_time<<endl;
    cout<<"xy vote time = "<<xy_vote_time<<endl;
    cout<<"z vote time = "<<z_vote_time<<endl;
    cout<<"cluster time = "<<cluster_time<<endl;
    cout<<"votesort time = "<<votesort_time<<endl;
  }
}


bool HelixHough::attemptClusterMerge(unsigned int zoomlevel, unsigned int MAX, unsigned int ca, unsigned int d, unsigned int r, unsigned int th, unsigned int zz0, unsigned int bin, unsigned int newbin, vector<unsigned char>& good_bins, unsigned int volume, float cluster_size_cut, float overlap_cut, vector<ParameterCluster>& clusters, unsigned int* bins_start, unsigned int* bins_end, vector<unsigned int>& map_clus, vector<unsigned char>& too_big, vector<unsigned int>& temp_merged, vector<unsigned int>& C)
{
  if((good_bins[newbin] == 1) && (map_clus[newbin] != 4294967295))
  {
    if(too_big[map_clus[newbin]] == 0)
    {
      unsigned int tempsize = clusters[map_clus[newbin]].hit_indexes.size();
      for(unsigned int ind=bins_start[bin];ind<=bins_end[bin];++ind)
      {
        clusters[map_clus[newbin]].hit_indexes.push_back((*(bins_vec[zoomlevel]))[ind].entry);
      }
      C.clear();
      C.assign(MAX+1, 0);
      in_place_counting_unique(clusters[map_clus[newbin]].hit_indexes, C);
      unsigned int size_diff = clusters[map_clus[newbin]].hit_indexes.size() - tempsize;
      unsigned int size_diff_2 = clusters[map_clus[newbin]].hit_indexes.size() - (1 + bins_end[bin] - bins_start[bin]);
      if( (((float)size_diff)/((float)(1 + bins_end[bin] - bins_start[bin])) < overlap_cut) || (((float)size_diff_2)/((float)(tempsize)) < overlap_cut) )
      {
        clusters[map_clus[newbin]].range.mergeRange(ca,d,r,th,zz0);
        ParameterCluster* cluster = &(clusters[map_clus[newbin]]);
        unsigned int cluster_volume = ((cluster->range.max_phi - cluster->range.min_phi)+1)*((cluster->range.max_d - cluster->range.min_d)+1)*((cluster->range.max_k - cluster->range.min_k)+1)*((cluster->range.max_dzdl - cluster->range.min_dzdl)+1)*((cluster->range.max_z0 - cluster->range.min_z0)+1);
        if( ((float)cluster_volume)/((float)volume) > cluster_size_cut )
        {
          too_big[map_clus[newbin]] = 1;
        }
        map_clus[bin] = map_clus[newbin];
        return true;
      }
      else
      {
        clusters[map_clus[newbin]].hit_indexes.resize(tempsize); 
      }
    }
  }
  return false;
}


void HelixHough::makeClusters(unsigned int zoomlevel, unsigned int MAX, unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0, unsigned int min_hits, vector<ParameterCluster>& clusters, bool& use_clusters, bool& is_super_bin)
{
  unsigned int volume = n_phi*n_d*n_k*n_dzdl*n_z0;
  float cluster_size_cut = 0.25;
  //float bin_size_cut = 0.75;
  float overlap_cut = 0.1;
  is_super_bin = false;
  
  vector<unsigned int> map_clus(volume, 4294967295);
  vector<unsigned char> good_bins(volume, 0);
  vector<unsigned char> too_big(volume, 0);
  
  for(unsigned int ca=0;ca<n_phi;++ca)
  {
    for(unsigned int d=0;d<n_d;++d)
    {
      for(unsigned int r=0;r<n_k;++r)
      {
        for(unsigned int th=0;th<n_dzdl;++th)
        {
          for(unsigned int zz0=0;zz0<n_z0;++zz0)
          {
            unsigned int bin = BinEntryPair5D::linearBin(n_d,n_k,n_dzdl,n_z0, ca,d,r,th,zz0);
            
            if(bins_end[bin] == 4294967295){continue;}
            if(( 1 + bins_end[bin] - bins_start[bin] ) >= min_hits )
            {
              if(check_layers == true)
              {
                unsigned int layer_mask = 0;
                for(unsigned int i=bins_start[bin];i<=bins_end[bin];++i)
                {
                  layer_mask = layer_mask | (1<<((*(hits_vec[zoomlevel]))[ (*(bins_vec[zoomlevel]))[i].entry ].layer));
                  unsigned int nlayers = __builtin_popcount(layer_mask); cout << "makeClusters_nlayers "<<nlayers<<endl;
                  if(nlayers>=req_layers){good_bins[bin] = 1;}
                }
              }
              else
              {
                good_bins[bin] = 1;
              }
            }
            else{continue;}
            
            if(good_bins[bin] == 0){continue;}
            
            if(ca > 0)
            {
              unsigned int newbin = bin - n_d*n_k*n_dzdl*n_z0;
              if(attemptClusterMerge(zoomlevel, MAX, ca, d, r, th, zz0, bin, newbin, good_bins, volume, cluster_size_cut, overlap_cut, clusters, (unsigned int*)(bins_start), (unsigned int*)(bins_end), map_clus, too_big, temp_merged_clus, C_clus) == true){clustersinzoomlevel_phi[zoomlevel]++;  continue;}
            }
            
            if(d > 0)
            {
              unsigned int newbin = bin - n_k*n_dzdl*n_z0;
              if(attemptClusterMerge(zoomlevel, MAX, ca, d, r, th, zz0, bin, newbin, good_bins, volume, cluster_size_cut, overlap_cut, clusters, (unsigned int*)(bins_start), (unsigned int*)(bins_end), map_clus, too_big, temp_merged_clus, C_clus) == true){clustersinzoomlevel_d[zoomlevel]++; continue;}
            }
            
            if(r > 0)
            {
              unsigned int newbin = bin - n_dzdl*n_z0;
              if(attemptClusterMerge(zoomlevel, MAX, ca, d, r, th, zz0, bin, newbin, good_bins, volume, cluster_size_cut, overlap_cut, clusters, (unsigned int*)(bins_start), (unsigned int*)(bins_end), map_clus, too_big, temp_merged_clus, C_clus) == true){clustersinzoomlevel_k[zoomlevel]++; continue;}
            }
            
            if(th > 0)
            {
              unsigned int newbin = bin - n_z0;
              if(attemptClusterMerge(zoomlevel, MAX, ca, d, r, th, zz0, bin, newbin, good_bins, volume, cluster_size_cut, overlap_cut, clusters, (unsigned int*)(bins_start), (unsigned int*)(bins_end), map_clus, too_big, temp_merged_clus, C_clus) == true){clustersinzoomlevel_dzdl[zoomlevel]++; continue;}
            }
            
            if(zz0 > 0)
            {
              unsigned int newbin = bin - 1;
              if(attemptClusterMerge(zoomlevel, MAX, ca, d, r, th, zz0, bin, newbin, good_bins, volume, cluster_size_cut, overlap_cut, clusters, (unsigned int*)(bins_start), (unsigned int*)(bins_end), map_clus, too_big, temp_merged_clus, C_clus) == true){clustersinzoomlevel_z0[zoomlevel]++; continue;}
            }
            if(num_clusters[zoomlevel] >= clusters.size())
            {
              clusters.push_back(ParameterCluster());
            }
            clusters[num_clusters[zoomlevel]].range = ParRange(ca,ca,d,d,r,r,th,th,zz0,zz0);
            map_clus[bin] = num_clusters[zoomlevel];
            clusters[num_clusters[zoomlevel]].hit_indexes.clear();
            for(unsigned int ind=bins_start[bin];ind<=bins_end[bin];++ind)
            {
              clusters[num_clusters[zoomlevel]].hit_indexes.push_back((*(bins_vec[zoomlevel]))[ind].entry);
            }
            
            num_clusters[zoomlevel] += 1;
          }
        }
      }
    }
  }
}


void HelixHough::findHelices(unsigned int min_hits, unsigned int max_hits, vector<SimpleTrack3D>& tracks, unsigned int maxtracks, unsigned int zoomlevel)
{
  unsigned int tracks_at_start = tracks.size();
  
  if((maxtracks != 0) && (tracks.size() >= max_tracks)){return;}
  
  hitsinzoomlevel[zoomlevel].push_back(hits_vec[zoomlevel]->size());




  timeval t1,t2;
  double time1=0.;
  double time2=0.;
  if(print_timings==true)
  {
    gettimeofday(&t1, NULL);
  }
  vote(zoomlevel);
  if(print_timings==true)
  {
    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
    vote_time += (time2 - time1);
  }
  
  unsigned int n_entries = bins_vec[zoomlevel]->size();
  if(n_entries == 0){return;}
  
  unsigned int n_phi = n_phi_bins[zoomlevel];
  unsigned int n_d = n_d_bins[zoomlevel];
  unsigned int n_k = n_k_bins[zoomlevel];
  unsigned int n_dzdl = n_dzdl_bins[zoomlevel];
  unsigned int n_z0 = n_z0_bins[zoomlevel];
  
  num_clusters[zoomlevel] = 0;
  bool use_clusters = true;
  bool is_super_bin = false;
  if(zoomlevel>=cluster_start_bin)
  {
    if(print_timings==true)
    {
      gettimeofday(&t1, NULL);
    }
    
    makeClusters(zoomlevel, hits_vec[zoomlevel]->size(), n_phi, n_d, n_k, n_dzdl, n_z0, min_hits, *(clusters_vec[zoomlevel]), use_clusters, is_super_bin);
    
    if(print_timings)
    {
      gettimeofday(&t2, NULL);
      time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
      time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
      cluster_time += (time2 - time1);
    }
  }
  else{use_clusters=false;}
  
  if(use_clusters == false)
  {
    unsigned int count = 0;
    hits_vec[zoomlevel+1]->clear();
    setRange((*(bins_vec[zoomlevel]))[count], zoomranges[zoomlevel], zoomranges[zoomlevel+1], n_phi, n_d, n_k, n_dzdl, n_z0);
    //scan over the bins in 5-D hough space
    while(count < n_entries)
    {
      if(remove_hits == true)
      {
        if((*hit_used)[(*(hits_vec[zoomlevel]))[ (*(bins_vec[zoomlevel]))[count].entry ].index] == false)
        {
          hits_vec[zoomlevel+1]->push_back( (*(hits_vec[zoomlevel]))[ (*(bins_vec[zoomlevel]))[count].entry ] );
        }
      }
      else
      {
        hits_vec[zoomlevel+1]->push_back( (*(hits_vec[zoomlevel]))[ (*(bins_vec[zoomlevel]))[count].entry ] );
      }
      
      
      
      count+=1;
      //we have collected all hits from this bin. now zoom again or find tracks with user routine
      if(  (count == n_entries) || ((*(bins_vec[zoomlevel]))[count].bin != (*(bins_vec[zoomlevel]))[count-1].bin)  )
      {
        if(hits_vec[zoomlevel+1]->size() >= min_hits)
        {
          if(breakRecursion(*(hits_vec[zoomlevel+1]), zoomranges[zoomlevel+1]) == true)
          {
          }
          else if( ((zoomlevel+1)==max_zoom) )
          {
            findTracks(*(hits_vec[zoomlevel+1]), tracks, zoomranges[zoomlevel+1]);
          }
          else if( ((zoomlevel+1) >= min_zoom) && ( (hits_vec[zoomlevel+1]->size() <= max_hits) ) )
          {
            findTracks(*(hits_vec[zoomlevel+1]), tracks, zoomranges[zoomlevel+1]);
          }
          else if( (hits_vec[zoomlevel+1]->size() <= max_hits_pairs) && ((zoomlevel+1) >= min_zoom) )
          {
            findHelicesByPairsBegin(min_hits, max_hits, tracks, maxtracks, zoomlevel+1);
          }
          else
          {
            findHelices(min_hits, max_hits, tracks, maxtracks, zoomlevel+1);
          }
          if(maxtracks != 0)
          {
            double phi_proportion = (zoomranges[zoomlevel].max_phi - zoomranges[zoomlevel].min_phi)/((zoomranges[0].max_phi - zoomranges[0].min_phi));
            double d_proportion = (zoomranges[zoomlevel].max_d - zoomranges[zoomlevel].min_d)/((zoomranges[0].max_d - zoomranges[0].min_d));
            double k_proportion = (zoomranges[zoomlevel].max_k - zoomranges[zoomlevel].min_k)/((zoomranges[0].max_k - zoomranges[0].min_k));
            unsigned int expected_tracks = (unsigned int)(fabs( ((double)maxtracks)*phi_proportion*d_proportion*k_proportion )) + 1;
            
            if((tracks.size() - tracks_at_start) > expected_tracks)
            {
              return;
            }
          }
        }
        if(count == n_entries){break;}
        hits_vec[zoomlevel+1]->clear();
        setRange((*(bins_vec[zoomlevel]))[count], zoomranges[zoomlevel], zoomranges[zoomlevel+1], n_phi, n_d, n_k, n_dzdl, n_z0);
      }
    }
  }
  else
  {
    if(clusters_vec[zoomlevel]->size() == 0){return;}
    // for each cluster, eiter perform the Hough again or break into user-defined routine
    for(unsigned int i=0, size=num_clusters[zoomlevel];i<size;++i)
    {
      hits_vec[zoomlevel+1]->clear();
      vector<unsigned int>::iterator index_iter;
      for(index_iter=(*(clusters_vec[zoomlevel]))[i].hit_indexes.begin();index_iter!=(*(clusters_vec[zoomlevel]))[i].hit_indexes.end();++index_iter)
      {
        if(remove_hits == true)
        {
          if( (*hit_used)[(*(hits_vec[zoomlevel]))[*index_iter].index] == false )
          {
            hits_vec[zoomlevel+1]->push_back((*(hits_vec[zoomlevel]))[*index_iter]);
          }
        }
        else
        {
          hits_vec[zoomlevel+1]->push_back((*(hits_vec[zoomlevel]))[*index_iter]);
        }
      }
      setClusterRange(zoomranges[zoomlevel], zoomranges[zoomlevel+1], (*(clusters_vec[zoomlevel]))[i].range, n_phi, n_d, n_k, n_dzdl, n_z0);
      if((breakRecursion(*(hits_vec[zoomlevel+1]), zoomranges[zoomlevel+1]) == true)){}
      else if((zoomlevel+1)==max_zoom)
      {
        findTracks(*(hits_vec[zoomlevel+1]), tracks, zoomranges[zoomlevel+1]);
      }
      else if( ((zoomlevel+1) >= min_zoom) && (hits_vec[zoomlevel+1]->size() <= max_hits) )
      {
        findTracks(*(hits_vec[zoomlevel+1]), tracks, zoomranges[zoomlevel+1]);
      }
      else if( (hits_vec[zoomlevel+1]->size() <= max_hits_pairs) && ((zoomlevel+1) >= min_zoom) )
      {
        findHelicesByPairsBegin(min_hits, max_hits, tracks, maxtracks, zoomlevel+1);
      }
      else
      {
        findHelices(min_hits, max_hits, tracks, maxtracks, zoomlevel+1);
      }
    }
  }
  
  
}


void HelixHough::findSeededHelices(vector<SimpleTrack3D>& seeds, vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits, vector<SimpleTrack3D>& tracks, unsigned int maxtracks)
{
  vote_time = 0.;
  xy_vote_time = 0.;
  z_vote_time = 0.;
  
  initEvent(hits, min_hits);
  initSeeding(seeds);
  
  hit_used->clear();
  index_mapping.clear();
  index_mapping.resize(hits.size(), 0);
  for(unsigned int i=0;i<hits.size();i++)
  {
    index_mapping[i] = hits[i].index;
    hits[i].index = i;
  }
  (*hit_used).assign(hits.size(), false);
  cout << "hits.size() "<<hits.size()<< " maxtracks "<<maxtracks<<" min_hits "<<min_hits<<endl;
  max_tracks = maxtracks;
  base_hits = &hits;
  (*(hits_vec[0])) = hits;cout << "hits size "<<hits.size()<<endl;
  (*(seeds_vec[0])) = seeds;
  current_range = top_range;
  zoomranges.clear();
  for(unsigned int z=0;z<=max_zoom;z++)
  {
    zoomranges.push_back(top_range);
  }
  vector<SimpleTrack3D> temp_tracks;
  
  if(separate_by_helicity==true)
  {
    helicity=true;
    findSeededHelices(min_hits, max_hits, temp_tracks, maxtracks, 0);
    helicity=false;
    findSeededHelices(min_hits, max_hits, temp_tracks, maxtracks, 0);
  }
  else{findSeededHelices(min_hits, max_hits, temp_tracks, maxtracks, 0);}
  
  for(unsigned int i=0;i<hits.size();i++)
  {
    hits[i].index = index_mapping[i];
  }
  
  finalize(temp_tracks, tracks);
  
  if(print_timings == true)
  {
    cout<<"vote time = "<<vote_time<<endl;
    cout<<"xy vote time = "<<xy_vote_time<<endl;
    cout<<"z vote time = "<<z_vote_time<<endl;
    cout<<"cluster time = "<<cluster_time<<endl;
    cout<<"votesort time = "<<votesort_time<<endl;
  }
}


class floatBin
{
public:
  floatBin(float l, float h) : low(l), high(h) {}
  ~floatBin(){}
  float low;
  float high;
  bool operator <(const floatBin& other) const
  {
    return ( high < other.low );
  }
};


//return which bin in bins that val belongs to, or -1 if it doesn't belong to any
static int seed_bin(vector<floatBin>& bins, float val)
{
  floatBin bin(val,val);
  if( ( bin < bins[0] ) || ( bins.back() < bin ) ){return -1;}
  pair<vector<floatBin>::iterator,vector<floatBin>::iterator> bounds;
  bounds = equal_range(bins.begin(), bins.end(), bin);
  return ( (int)(bounds.first - bins.begin()) );
}


void HelixHough::setRangeFromSeed(HelixRange& range, SimpleTrack3D& seed)
{
  range.min_phi = seed.phi - 0.001;
  range.max_phi = seed.phi + 0.001;
  if(range.min_phi < 0.){range.min_phi = 0.;}
  if(range.max_phi > 2.*M_PI){range.max_phi = 2.*M_PI;}
  range.min_d = seed.d - 0.001;
  range.max_d = seed.d + 0.001;
  range.min_k = seed.kappa - 0.00001;
  range.max_k = seed.kappa + 0.00001;
  if(range.min_k < 0.){range.min_k = 0.;}
  
  range.min_k = range.min_k * range.min_k;
  range.max_k = range.max_k * range.max_k;
  
  range.min_dzdl = seed.dzdl - 0.05;
  range.max_dzdl = seed.dzdl + 0.05;
  range.min_z0 = seed.z0 - 0.01;
  range.max_z0 = seed.z0 + 0.01;
}


void HelixHough::findSeededHelices(unsigned int min_hits, unsigned int max_hits, vector<SimpleTrack3D>& tracks, unsigned int maxtracks, unsigned int zoomlevel)
{
  unsigned int tracks_at_start = tracks.size();
  
  if((maxtracks != 0) && (tracks.size() >= max_tracks)){return;}
  bool debug = false;
   
  timeval t1,t2,t3,t4;
  double time1=0.;
  double time2=0.;
  double time3=0.;
  double time4=0.;
  gettimeofday(&t1, NULL);
  vote(zoomlevel);
//moved here....
  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
  vote_time += (time2 - time1);
  
  unsigned int n_phi = n_phi_bins[zoomlevel];
  unsigned int n_d = n_d_bins[zoomlevel];
  unsigned int n_k = n_k_bins[zoomlevel];
  unsigned int n_dzdl = n_dzdl_bins[zoomlevel];
  unsigned int n_z0 = n_z0_bins[zoomlevel];
  
  float low=0.;
  float high=0.;
  float delta=0.;
  if (debug) cout << "hits_vec["<<zoomlevel<<"] "<<hits_vec[zoomlevel]->size()<<"seeds_vec["<<zoomlevel<<"]->size()"<<seeds_vec[zoomlevel]->size()<<endl;
  vector<floatBin> phibins;
  low = zoomranges[zoomlevel].min_phi;
  delta = (zoomranges[zoomlevel].max_phi - zoomranges[zoomlevel].min_phi)/((float)(n_phi));
  high = low+delta;
  for(unsigned int i=0;i<n_phi;++i)
  {
    phibins.push_back(floatBin(low,high));
    low+=delta;
    high+=delta;
  }
  
  vector<floatBin> dbins;
  low = zoomranges[zoomlevel].min_d;
  delta = (zoomranges[zoomlevel].max_d - zoomranges[zoomlevel].min_d)/((float)(n_d));
  high = low+delta;
  for(unsigned int i=0;i<n_d;++i)
  {
    dbins.push_back(floatBin(low,high));
    low+=delta;
    high+=delta;
  }
  
  vector<floatBin> kbins;
  low = zoomranges[zoomlevel].min_k;
  delta = (zoomranges[zoomlevel].max_k - zoomranges[zoomlevel].min_k)/((float)(n_k));
  high = low+delta;
  for(unsigned int i=0;i<n_k;++i)
  {
    kbins.push_back(floatBin(low,high));
    low+=delta;
    high+=delta;
  }
  
  vector<floatBin> z0bins;
  low = zoomranges[zoomlevel].min_z0;
  delta = (zoomranges[zoomlevel].max_z0 - zoomranges[zoomlevel].min_z0)/((float)(n_z0));
  high = low+delta;
  for(unsigned int i=0;i<n_z0;++i)
  {
    z0bins.push_back(floatBin(low,high));
    low+=delta;
    high+=delta;
  }
  
  vector<floatBin> dzdlbins;
  low = zoomranges[zoomlevel].min_dzdl;
  delta = (zoomranges[zoomlevel].max_dzdl - zoomranges[zoomlevel].min_dzdl)/((float)(n_dzdl));
  high = low+delta;
  for(unsigned int i=0;i<n_dzdl;++i)
  {
    dzdlbins.push_back(floatBin(low,high));
    low+=delta;
    high+=delta;
  }
  
  if (debug) cout << "seeds_vec["<<zoomlevel<<"] "<<seeds_vec[zoomlevel]->size()<<endl;
  // voting for the seeds
  for(unsigned int i=0;i<seeds_vec[zoomlevel]->size();++i)
  {
    float d = (*(seeds_vec[zoomlevel]))[i].d;
    float phi = (*(seeds_vec[zoomlevel]))[i].phi;
    float kappa = pow((*(seeds_vec[zoomlevel]))[i].kappa, 1.);
    float z0 = (*(seeds_vec[zoomlevel]))[i].z0;
    float dzdl = (*(seeds_vec[zoomlevel]))[i].dzdl;
    
    int tempbin = 0;
    
    unsigned int phi_bin = 0;
    tempbin = seed_bin(phibins, phi);
    if(tempbin < 0){continue;}
    else{phi_bin = (unsigned int)(tempbin);}
    
    unsigned int d_bin = 0;
    tempbin = seed_bin(dbins, d);
    if(tempbin < 0){continue;}
    else{d_bin = (unsigned int)(tempbin);}
    
    unsigned int kappa_bin = 0;
    tempbin = seed_bin(kbins, kappa);
    if(tempbin < 0){continue;}
    else{kappa_bin = (unsigned int)(tempbin);}
    
    unsigned int z0_bin = 0;
    tempbin = seed_bin(z0bins, z0);
    if(tempbin < 0){continue;}
    else{z0_bin = (unsigned int)(tempbin);}
    
    unsigned int dzdl_bin = 0;
    tempbin = seed_bin(dzdlbins, dzdl);
    if(tempbin < 0){continue;}
    else{dzdl_bin = (unsigned int)(tempbin);}
    
    unsigned int bin = BinEntryPair5D::linearBin(n_d, n_k, n_dzdl, n_z0, phi_bin, d_bin, kappa_bin, dzdl_bin, z0_bin);
    
    bins_vec[zoomlevel]->push_back(BinEntryPair5D(bin, i, true));
  }
  //int seedshere = 0;
  //for(unsigned int i=0;i<bins_vec[zoomlevel]->size();i++)
  //if ((*(bins_vec[zoomlevel]))[i].is_seed) seedshere++;
  //cout <<"seedshere 1 = "<<seedshere;
  if (debug) cout << "bins_vec["<<zoomlevel<<"] "<<bins_vec[zoomlevel]->size()<<endl;
  gettimeofday(&t3, NULL);
  sort(bins_vec[zoomlevel]->begin(), bins_vec[zoomlevel]->end());
  gettimeofday(&t4, NULL);
  time3 = ((double)(t3.tv_sec) + (double)(t3.tv_usec)/1000000.);
  time4 = ((double)(t4.tv_sec) + (double)(t4.tv_usec)/1000000.);
  votesort_time += (time4 - time3);
  

  unsigned int n_entries = bins_vec[zoomlevel]->size();
  if(n_entries == 0){return;}

  //seedshere = 0;
  //for(unsigned int i=0;i<bins_vec[zoomlevel]->size();i++)
  //if ((*(bins_vec[zoomlevel]))[i].is_seed) seedshere++;
  //cout <<"seedshere 2 = "<<seedshere;


  unsigned int count = 0;
  hits_vec[zoomlevel+1]->clear();
  seeds_vec[zoomlevel+1]->clear();
  setRange((*(bins_vec[zoomlevel]))[count], zoomranges[zoomlevel], zoomranges[zoomlevel+1], n_phi, n_d, n_k, n_dzdl, n_z0);
  //scan over the bins in 5-D hough space

  //seedshere = 0;
  //for(unsigned int i=0;i<bins_vec[zoomlevel]->size();i++)
  //if ((*(bins_vec[zoomlevel]))[i].is_seed) seedshere++;
  //cout <<"seedshere 3 = "<<seedshere;


  while(count < n_entries)
  {
    if( (*(bins_vec[zoomlevel]))[count].is_seed == false)
    {//cout <<"is_seed"<<endl;
      if(remove_hits == true)
      {
        if((*hit_used)[(*(hits_vec[zoomlevel]))[ (*(bins_vec[zoomlevel]))[count].entry ].index] == false)
        {
          hits_vec[zoomlevel+1]->push_back( (*(hits_vec[zoomlevel]))[ (*(bins_vec[zoomlevel]))[count].entry ] );
        }
      }
      else
      {//cout <<"one new!" <<endl;
        hits_vec[zoomlevel+1]->push_back( (*(hits_vec[zoomlevel]))[ (*(bins_vec[zoomlevel]))[count].entry ] );
      }
    }
    else
    {
      seeds_vec[zoomlevel+1]->push_back( (*(seeds_vec[zoomlevel]))[ (*(bins_vec[zoomlevel]))[count].entry ] );
    }
    //cout << "3!"<<count<<endl;

    
    count++;
    //we have collected all hits from this bin. now zoom again or find tracks with user routine
    if(  (count == n_entries) || ((*(bins_vec[zoomlevel]))[count].bin != (*(bins_vec[zoomlevel]))[count-1].bin)  )
    { //cout << "4! "<<hits_vec[zoomlevel+1]->size()<<endl;
      if( (hits_vec[zoomlevel+1]->size() >= min_hits) && (seeds_vec[zoomlevel+1]->size() != 0) )
	{//cout <<"BR??"<<endl;
        if(breakRecursion(*(hits_vec[zoomlevel+1]), zoomranges[zoomlevel+1]) == true)
	  {//cout <<"BR zoomlevel"<<zoomlevel<<endl;
          findSeededTracks(*(seeds_vec[zoomlevel+1]), *(hits_vec[zoomlevel+1]), tracks, zoomranges[zoomlevel+1]);
        }
        else if( ((zoomlevel+1) >= min_zoom) && ( (hits_vec[zoomlevel+1]->size() <= max_hits) || (zoomlevel+1)==max_zoom ) )
	  {if (debug) cout <<"maxhits "<<max_hits<<" hits_vec[zoomlevel+1]->size()"<<hits_vec[zoomlevel+1]->size()<<" zoomlevel "<<zoomlevel<<" (zoomlevel+1)==max_zoom "<<(((zoomlevel+1)==max_zoom)?"true":"false")<<endl;
          findSeededTracks(*(seeds_vec[zoomlevel+1]), *(hits_vec[zoomlevel+1]), tracks, zoomranges[zoomlevel+1]);
        }
        else
	  { if (debug) cout <<"inside!"<<endl;
          findSeededHelices(min_hits, max_hits, tracks, maxtracks, zoomlevel+1);
        }
        if(maxtracks != 0)
        {
          double phi_proportion = (zoomranges[zoomlevel].max_phi - zoomranges[zoomlevel].min_phi)/((zoomranges[0].max_phi - zoomranges[0].min_phi));
          double d_proportion = (zoomranges[zoomlevel].max_d - zoomranges[zoomlevel].min_d)/((zoomranges[0].max_d - zoomranges[0].min_d));
          double k_proportion = (zoomranges[zoomlevel].max_k - zoomranges[zoomlevel].min_k)/((zoomranges[0].max_k - zoomranges[0].min_k));
          unsigned int expected_tracks = (unsigned int)(fabs( ((double)maxtracks)*phi_proportion*d_proportion*k_proportion )) + 1;
          
          if((tracks.size() - tracks_at_start) > expected_tracks)
          {cout <<"WHAT?!"<<endl;
            return; 
          }
        }
      }
      if(count == n_entries){break;}
      hits_vec[zoomlevel+1]->clear();
      seeds_vec[zoomlevel+1]->clear();
      setRange((*(bins_vec[zoomlevel]))[count], zoomranges[zoomlevel], zoomranges[zoomlevel+1], n_phi, n_d, n_k, n_dzdl, n_z0);
    }
    //std::cout<<"out of loop"<<std::endl;

  }
}



