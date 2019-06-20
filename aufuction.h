namespace aufunc
{
	template<class Type>
	inline bool bNormal(const std::pair<Type,Type>* quad)
	{
		Type k[4],Langle[4],angle[4];
		const Type PI=3.1415926;
		Type dx,dy,d;
		//求每一条线的角度
		for(size_t i=0;i<3;++i)
		{
			dx=quad[i+1].first-quad[i].first;
			dy=quad[i+1].second-quad[i].second;
			d=std::sqrt(dx*dx+dy*dy);
			if(d<1e-6) return false;
			Langle[i]=std::asin(dy/d);
			if(dx<0) Langle[i]=PI-Langle[i];	
		}
		dx=quad[0].first-quad[3].first;
		dy=quad[0].second-quad[3].second;
		d=std::sqrt(dx*dx+dy*dy);
		if(d<1e-6) return false;
		Langle[3]=std::asin(dy/d);
		if(dx<0) Langle[3]=PI-Langle[3];
		//求每一个角的角度
		for(size_t i=1;i<=3;i++)
		{
			Type tmp=Langle[i-1]-Langle[i];
			if(tmp>=PI) 
			{
				tmp-=2*PI;
			}
			if(tmp<-PI) 
			{
				tmp+=2*PI;
			}
			angle[i]=PI+tmp;
		}
		Type tmp=Langle[3]-Langle[0];
			if(tmp>=PI) 
			{
				tmp-=2*PI;
			}
			if(tmp<-PI) 
			{
				tmp+=2*PI;
				Langle[0]-=2*PI;
			}
		angle[0]=PI+tmp;
		for (size_t i=0;i<4;i++)
			angle[i]=angle[i]*180/PI;
		tmp=angle[0]+angle[1]+angle[2]+angle[3];
		if(std::abs(tmp-720)<1) return false;
		if(abs(tmp-360)>1)
		{
			for (size_t i=0;i<4;i++)
				if((360-angle[i])>=180||(360-angle[i])<1e-6) return false;
		}
		else
		{
			for (size_t i=0;i<4;i++)
				if(angle[i]>=180||angle[i]<1e-6) return false;
		}
		return true;
	}
}

namespace aufunc
{
	template<class Type>
	inline bool bInArea(const std::pair<Type,Type>* quad,const std::pair<Type,Type>point)
	{
		if(!bNormal(quad)) return false;
		std::complex<Type> vec1[4],vec2[4];
		for(size_t i=0;i<3;++i)
		{
			vec1[i]=std::complex<Type>(quad[i+1].first-quad[i].first,quad[i+1].second-quad[i].second);
			vec2[i]=std::complex<Type>(point.first-quad[i].first,point.second-quad[i].second);
		}
		vec1[3]=std::complex<Type>(quad[0].first-quad[3].first,quad[0].second-quad[3].second);
		vec2[3]=std::complex<Type>(point.first-quad[3].first,point.second-quad[3].second);
		int flag;
		Type angle;
		if(vec2[0]==Type(0)) return true;
		angle=std::arg(vec1[0]/vec2[0]);
		if(angle==0) flag=0;
		else if(angle>0) flag=1;
		else flag=-1;
		for(size_t i=1;i<4;++i)
		{
			if(vec2[i]==Type(0)) return true;
			angle=std::arg(vec1[i]/vec2[i]);
			if(angle>0) 
			{
				if(flag>0) continue;
				else if(flag==0)
				{
					flag=1;
					continue;
				}
				else return false;
			}
			else if(angle<0)
			{
				if(flag<0) continue;
				else if(flag==0)
				{
					flag=-1;
					continue;
				}
				else return false;
			}
		}
		return true;
	}
}

namespace aufunc
{
	template<class Type>
	inline bool bInAreaT(const std::pair<Type,Type>* tri,const std::pair<Type,Type>point)
	{
		std::complex<Type> vec1[4],vec2[4];
		for(size_t i=0;i<2;++i)
		{
			vec1[i]=std::complex<Type>(tri[i+1].first-tri[i].first,tri[i+1].second-tri[i].second);
			vec2[i]=std::complex<Type>(point.first-tri[i].first,point.second-tri[i].second);
		}
		vec1[2]=std::complex<Type>(tri[0].first-tri[2].first,tri[0].second-tri[2].second);
		vec2[2]=std::complex<Type>(point.first-tri[2].first,point.second-tri[2].second);
		int flag;
		Type angle;
		if(vec2[0]==Type(0)) return true;
		angle=std::arg(vec1[0]/vec2[0]);
		if(angle==0) flag=0;
		else if(angle>0) flag=1;
		else flag=-1;
		for(size_t i=1;i<3;++i)
		{
			if(vec2[i]==Type(0)) return true;
			angle=std::arg(vec1[i]/vec2[i]);
			if(angle>0) 
			{
				if(flag>0) continue;
				else if(flag==0)
				{
					flag=1;
					continue;
				}
				else return false;
			}
			else if(angle<0)
			{
				if(flag<0) continue;
				else if(flag==0)
				{
					flag=-1;
					continue;
				}
				else return false;
			}
		}
		return true;
	}
}

namespace aufunc
{
	template<class Type>
	inline Type conj(const Type &orig)
	{
		return orig;
	}

	template<class Type>
	inline std::complex<Type> conj(const std::complex<Type> &orig)
	{
		return std::conj(orig);
	}
}

namespace aufunc
{
	template<class Type>
	inline Type abs(const Type &orig)
	{
		return std::abs(orig);
	};

	template<class Type>
	inline long double abs(const std::complex<Type> &orig)
	{
		long double tmp1=std::real(orig),tmp2=std::imag(orig);
		return std::pow(tmp1*tmp1+tmp2*tmp2,long double(0.5));
	}
}

namespace aufunc
{
	template<class Type>
	inline Type real(const Type &orig)
	{
		return orig;
	};

	template<class Type>
	inline Type real(const std::complex<Type> &orig)
	{
		return orig.real();
	}

	template<class Type>
	inline Type imag(const Type &orig)
	{
		return Type(0);
	};

	template<class Type>
	inline Type imag(const std::complex<Type> &orig)
	{
		return orig.imag();
	}
}

namespace aufunc
{
	template<class Type>
	Type rand(Type v1,Type v2)
	{
		if(v1==v2) return v1;
		if(v1>v2) std::swap(v1,v2);
		std::srand(time(0));
		std::rand();
		const int value=std::rand();
		const int ys=value%(v2-v1);
		const int cs=(value-ys)/(v2-v1);
		if(ys!=0)
			return ys+v1;
		else
		{
			if(cs%2==0) return v1;
			else return v2;
		}
	}

	template<class Type>
	Type randf(Type v1,Type v2)
	{
		std::srand(time(0));
		std::rand();
		return (v2-v1)/Type(RAND_MAX)*Type(std::rand())+v1;
	}
}