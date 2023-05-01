//------------------------------------------------------------------------------
//  Copyright 2007-2012 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "acd2d.h"
#include "acd2d_util.h"
#include "acd2d_cut.h"
#include "acd2d_dir.h"

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

namespace acd2d
{
	///////////////////////////////////////////////////////////////////////////
	
	cd_2d::cd_2d(bool save_diagonal)
	{
		store_diagoanls=save_diagonal;
		alpha=0; beta=1;
	}
	
	cd_2d::~cd_2d()
	{
		destroy();
	}
	
	///////////////////////////////////////////////////////////////////////////
	// polygon functions
	
	void cd_2d::addPolygon(const cd_polygon& poly)
	{
		if(poly.valid())
		{
			cd_polygon mypoly;
			todo_list.push_back(mypoly);
			todo_list.back().copy(buf_, poly);
			todo_list.back().buildDependency();
		}
		else
			cerr<<"! Error: acd_2d::addPolygon: Not a valid polygon"<<endl;
	}
	
	void cd_2d::destroy()
	{
		todo_list.clear();
		done_list.clear();
	}
	
	///////////////////////////////////////////////////////////////////////////
	void cd_2d::decomposeAll(double d, IConcavityMeasure * measure)
	{	const int iter_limit = 100;
		int iter = 0;
		if( d<1e-20 ) d=1e-20;
		do {
			try {
				decompose(d,measure);
				iter++;
			} catch(const std::exception& e) {
				throw std::runtime_error(e.what());
			}			
		}
		while(!todo_list.empty() and iter < iter_limit);
	}
    void cd_2d::maybe_decomposeAll(double d, IConcavityMeasure * measure)
    {	const int iter_limit = 100;
      int iter = 0;
      if( d<1e-20 ) d=1e-20;
      do {
          maybe_decompose(d, measure);
          iter++;
      }
      while(!todo_list.empty() and iter < iter_limit);
    }

	void cd_2d::decompose(double d, IConcavityMeasure * measure)
	{
		list<cd_polygon> ps;
		ps.swap(todo_list);
		list<cd_polygon>::iterator ips=ps.begin();
		m_measure=measure;
		if( m_measure==NULL ) {
			cerr<<"! ERROR: cd_2d::decompose: measure si NULL"<<endl;
			return;
		}
		if( d<1e-20 ) d=1e-20;
	
		for(;ips!=ps.end();ips++){
			cd_polygon& polys=*ips;
			try {
				decompose(d,polys);
			} catch (exception e) {
				throw std::runtime_error(e.what());
			}
			
		}
	}

    void cd_2d::maybe_decompose(double d, IConcavityMeasure * measure)
    {
      list<cd_polygon> ps;
      ps.swap(todo_list);
      list<cd_polygon>::iterator ips=ps.begin();
      m_measure=measure;
      if( m_measure==nullptr ) {
        cerr<<"! ERROR: cd_2d::decompose: measure si NULL"<<endl;
        return;
      }
      if( d<1e-20 ) d=1e-20;

      for(;ips!=ps.end();ips++){
        cd_polygon& polys=*ips;
        maybe_decompose(d,polys);
      }
    }
	
	void cd_2d::decompose(double d, cd_polygon& polys) {
		//if there are inner polys, random pick one and find the cut
		cd_poly poly=polys.next();
		if( poly.getType()==cd_poly::PIN ) // hole
			decompose_IN(d,polys,poly);
		else {
			try {
				decompose_OUT(d,polys,findOutMost(polys));
			} catch (exception e) {
				throw std::runtime_error(e.what());
			}
		}
	}

    void cd_2d::maybe_decompose(double d, cd_polygon& polys) {
      //if there are inner polys, random pick one and find the cut
      cd_poly poly = polys.next();
      if( poly.getType()==cd_poly::PIN ) // hole
        decompose_IN(d,polys,poly);
      else {
        maybe_decompose_OUT(d,polys,findOutMost(polys));
      }
    }
	
	void cd_2d::decompose_OUT(double d, cd_polygon& polys, cd_poly& poly)
	{
		cd_line cut_l; //cut line
	
		//check if we need to cut it.
		cd_vertex * r=poly.findCW(m_measure).first;
	
		if( r==NULL ){
			done_list.push_back(polys);
			return;
		}
	
		//smaller than tolerance
		if( r->getConcavity()<d ){
			done_list.push_back(polys);
			return;
		}
	
		find_a_good_cutline(cut_l,r,alpha,beta);
	
		//cut into two polys
		pair<cd_polygon,cd_polygon> sub_polys;
		try{
			cd_diagonal dia=cutPolys(sub_polys,polys.front(),cut_l, buf_);
			if(store_diagoanls) dia_list.push_back(dia);
		} catch (exception e) {
			throw std::runtime_error(e.what());
		}
	
		//add into to do
		todo_list.push_back(sub_polys.first);
		todo_list.push_back(sub_polys.second);
	}


    void cd_2d::maybe_decompose_OUT(double d, cd_polygon& polys, cd_poly& poly)
    {
      cd_line cut_l; //cut line

      //check if we need to cut it.
      cd_vertex * r=poly.findCW(m_measure).first;

      if( r==nullptr ){
        done_list.push_back(polys);
        return;
      }

      //smaller than tolerance
      if( r->getConcavity()<d ){
        done_list.push_back(polys);
        return;
      }

      find_a_good_cutline(cut_l,r,alpha,beta);

      //cut into two polys
      pair<cd_polygon,cd_polygon> sub_polys;
      try{
        cd_diagonal dia=cutPolys(sub_polys,polys.front(),cut_l, buf_);
        if(store_diagoanls) dia_list.push_back(dia);
        //add into to do
        todo_list.push_back(sub_polys.first);
        todo_list.push_back(sub_polys.second);
      } catch (exception e) {
        return;
      }
    }

	void cd_2d::decompose_IN(double d, cd_polygon& polys, cd_poly& poly)
	{
		//find the out most boundary
		cd_poly& out=findOutMost(polys);
	
		//find concavity witness
		cd_vertex * r=poly.getCW().first;
	
		//find which cw is better
		cd_line cut_l;
		find_a_good_cutline_for_hole(cut_l,r,out);
		cd_diagonal dia=mergeHole(out,poly,cut_l, buf_);
		todo_list.push_back(polys);
	
		//store cut line
		if(store_diagoanls) dia_list.push_back(dia);
	}

}//namespace acd2d