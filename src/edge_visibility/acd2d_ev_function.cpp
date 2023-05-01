//------------------------------------------------------------------------------
//  Copyright 2007-2012 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "acd2d_ev_function.h"
#include "earcut.hpp"
#include <vector>
#include <array>

using std::vector;

namespace acd2d {
	
	typedef pair<int,ev_triangle *> Edge; //a pair of pt id and tri id
	typedef vector<Edge> EV;   //vector of Edge
	
	inline void add(EV& ev, int pid, ev_triangle * t)
	{
		Edge dte(pid,t);
		ev.push_back(dte);
	}
	
	inline ev_triangle * find(EV& ev, int pid) //return tri
	{
		typedef EV::iterator it;
		for(it i=ev.begin();i!=ev.end();i++)
			if( i->first==pid ) return i->second;
			return NULL;
	}
	
	//returns the triangle contains (0,N-1) and build the dual of given triangles
	inline ev_triangle * dualT(int * t, int tsize, int psize,  ev_tri_buffer& buf )
	{
		int last=psize-1;
	
		ev_triangle * start_T=NULL; //the triangle with edge (0,psize-1)
		ev_triangle * first_T=NULL;
	
		vector<EV> hash(psize,EV());
		
		for( int it=0;it<tsize;it++ ){ //for each triangle
			ev_triangle * cur_tri = buf.getNew();
			if( it==0 ) first_T=cur_tri;
			for(int ie=0;ie<3;ie++ ){ //for each edge in this triangle
				int vid1=t[it*3+ie];       //edge start
				int vid2=t[it*3+(ie+1)%3]; //edge end
				cur_tri->v[ie]=vid1;
				if( (vid1==0&&vid2==last) || (vid1==last&&vid2==0) ){ 
					start_T=cur_tri; 
				}
	
				ev_triangle * nei_tri=find(hash[vid1],vid2); //neighboring tri
				if( nei_tri==NULL ){ //not found
					add(hash[vid2],vid1,cur_tri); continue; 
				}
				//found
				{for(int i=0;i<3;i++) if(cur_tri->t[i]==NULL){cur_tri->t[i]=nei_tri; break;}}
				{for(int i=0;i<3;i++) if(nei_tri->t[i]==NULL){nei_tri->t[i]=cur_tri; break;}}
			}
		}//end it
		return start_T;
	}
	
	//return the triangle contains s and e
	ev_triangle * triangulate( ev_vertex * pts, int polysize,  ev_tri_buffer& buf)
	{
		//prepare for triangulation
		vector<vector<std::array<double, 2>>> v = {vector<std::array<double, 2>>(polysize, {0,0})};
		
		//copy vertices
		for(int i=0;i<polysize;i++){
			const Point2d& pos=pts[i].v->getPos();
			std::copy(pos.get(), pos.get() + 2, std::begin(v.front().at(i)));
		}//end for
		
		auto t = mapbox::earcut<int>(v);

        if (t.size() > 0) {
          return dualT(t.data(), t.size() / 3, polysize, buf);
        }
        return nullptr;
	}

	
} //namespace acd2d
