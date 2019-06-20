//该文件为Cmatrix的实现文件

//1.下标检查
// a.operator()仅进行非0检查,越界检查交给vector
// b.构造函数及其它内存操作函数不进行下标检查
// c.其它所有函数(包括拷贝构造函数)均进行下标非零和越界检查
// d.所有下标检查失败均抛出异常
//2.空矩阵检查
// a.仅对数值类的矩阵运算进行空矩阵检查
// b.对能返回计算成败信息的,直接返回错误代码,其它抛出异常
//3.目标矩阵和源矩阵相同性检查
// a.对拷贝构造进行检查抛出异常
// b.对直接赋值,交换进行检查以提高效率,但不抛出异常


template<class Type>
long double Cmatrix<Type>::precision=1e-6;
//构造函数*******************************************

template <class Type>
inline Cmatrix<Type>::Cmatrix(size_t row,size_t line, Type elem)
{
	row_count=row;
	line_count=line;
	matrix.resize(line_count);
    for (size_type i=0;i<line_count;matrix[i++].resize(row_count,elem));

};

template <class Type>
inline Cmatrix<Type>::Cmatrix(Method flag,size_type row)
{
	row_count=row;
	line_count=row;
	matrix.resize(line_count);
    for (size_type i=0;i<line_count;matrix[i++].resize(row_count,Type()));
	if(flag==UNIT)
	for (size_t i=1;i<=row;i++)
		for (size_t j=1;j<=row;j++)
			if(i==j) matrix[j-1][i-1]=1;
			else matrix[j-1][i-1]=0;
}

template <class Type>
inline Cmatrix<Type>::Cmatrix(Method flag,Type v1,Type v2,size_t row,size_t line)
{
	row_count=row;
	line_count=line;
	matrix.resize(line_count);
    for (size_type i=0;i<line_count;matrix[i++].resize(row_count,Type()));
	rand(v1,v2,flag);
}

template<class Type>
template<class U>
inline Cmatrix<Type>::Cmatrix(const U* ptr,size_type num,size_type row,size_type line,Type elem,Method flag)
{
	row_count=row;
	line_count=line;
	matrix.resize(line_count);
    for (size_type i=0;i<line_count;matrix[i++].resize(row_count,elem));
	assign(ptr,num,flag);
};


template <class Type>
template <class	U>
inline Cmatrix<Type>::Cmatrix(const Cmatrix<U>& orig)
{
	if(this==&orig) throw bad_oper("源矩阵与目标矩阵为同一矩阵!");
	*this=orig;
};

template<class Type>
template<class U>
inline Cmatrix<Type>::Cmatrix(const Cmatrix<U> &orig, size_t beg, size_t end,Method flag)
{
	if(beg==0||end==0) throw bad_oper("构造函数出错:行列标不能为0!");
	if(this==&orig) throw bad_oper("构造函数出错:源矩阵与目标矩阵为同一矩阵!");
	if(flag==LINE)
	{
		if(end>orig.line()||beg>orig.line()) throw bad_oper("构造函数出错:列标越界!");
		if(beg>end) std::swap(beg,end);
		row_count=orig.row();
		line_count=end-beg+1;
		matrix.resize(line_count);
		for (size_type i=0;i<line_count;matrix[i++].resize(row_count));
		for(size_t i=1;i<=row_count;i++)
			for(size_t j=beg;j<=end;j++)
				(*this)(i,j-beg+1)=orig(i,j);
	}
	else
	{
		if(end>orig.row()||beg>orig.row()) throw bad_oper("构造函数出错:行标越界!");
		if(beg>end) std::swap(beg,end);
		row_count=end-beg+1;
		line_count=orig.line();
		matrix.resize(line_count);
		for (size_type i=0;i<line_count;matrix[i++].resize(row_count));
		for(size_t i=beg;i<=end;i++)
			for(size_t j=1;j<=line_count;j++)
				(*this)(i-beg+1,j)=orig(i,j);
	}
}

template<class Type>
template<class U>
inline Cmatrix<Type>::Cmatrix(const Cmatrix<U> &orig, size_t row1, size_t line1,size_t row2,size_t line2)
{
	if(row1==0||line1==0||row2==0||line2==0) throw bad_oper("构造函数出错:行列标不能为0!");
	if(this==&orig) throw bad_oper("构造函数出错:源矩阵与目标矩阵为同一矩阵!");
	if(line1>orig.line()||line2>orig.line()) throw bad_oper("构造函数出错:列标越界!");
	if(row1>orig.row()||row2>orig.row()) throw bad_oper("构造函数出错:行标越界!");
	if(row1>row2) std::swap(row1,row2);
	if(line1>line2) std::swap(line1,line2);
	row_count=row2-row1+1;
	line_count=line2-line1+1;
	matrix.resize(line_count);
	for (size_type i=0;i<line_count;matrix[i++].resize(row_count));
	for(size_t i=row1;i<=row2;i++)
		for(size_t j=line1;j<=line2;j++)
			(*this)(i-row1+1,j-line1+1)=orig(i,j);
}

template<class Type>
template<class InputIterator>
inline Cmatrix<Type>::Cmatrix(InputIterator beg,InputIterator end,size_type row,size_type line,const Type &elem,Method flag)
{
	row_count=row;
	line_count=line;
	matrix.resize(line_count);
    for (size_type i=0;i<line_count;matrix[i++].resize(row_count,elem));
	if (flag==LINE)
	{
		for(size_type i=0;i<line;++i)
			for (size_type j=0;j<row;++j)
				if(beg!=end) matrix[i][j]=*beg++;
	}
	else
	{

		for(size_type i=0;i<row;++i)
			for (size_type j=0;j<line;++j)
				if(beg!=end) matrix[j][i]=*beg++;
	}
};



//非变动型矩阵查询************************************************
template<class Type>
inline bool Cmatrix<Type>::IsEmpty() const
{return line_count==0||row_count==0;}

template<class Type>
inline const Type &Cmatrix<Type>::operator ()(size_type row,size_type line) const
{
	if(row==0)throw bad_oper("operator()出错:行标不能为0!");
	if(line==0)throw bad_oper("operator()出错:列标不能为0!");
	return matrix[line-1][row-1];
	
};

template <class Type>
inline size_t Cmatrix<Type>::line() const
{return line_count;};

template <class Type>
inline size_t Cmatrix<Type>::row() const
{return row_count;};

template <class Type>
inline std::pair<size_t,size_t> Cmatrix<Type>::size() const
{return std::make_pair(row_count,line_count);};

template <class Type>
inline Cmatrix<Type>  Cmatrix<Type>::getrow(size_t start,size_t end) const
{
	if(start==0||end==0) throw bad_oper("getrow出错:行标不能为0!");
	if(start>end) std::swap(start,end);
	if (end>row_count) throw bad_oper("getrow出错:行标越界!");
	Cmatrix<Type> tmp(end-start+1,line_count);
	for(size_t i=1;i<=tmp.row();i++)
		for (size_t j=1;j<=line_count;j++)
			tmp(i,j)=(*this)(start+i-1,j);
	return tmp;
};

template <class Type>
inline Cmatrix<Type>  Cmatrix<Type>::getline(size_t start,size_t end) const
{
	if(start==0||end==0) throw bad_oper("getline出错:列标不能为0!");
	if(start>end) std::swap(start,end);
	if (end>line_count) throw bad_oper("getline出错:列标越界!");
	Cmatrix<Type> tmp(row_count,end-start+1);
	for(size_t i=1;i<=row_count;i++)
		for (size_t j=1;j<=tmp.line();j++)
			tmp(i,j)=(*this)(i,start+j-1);
	return tmp;
};

template <class Type>
inline Cmatrix<Type>  Cmatrix<Type>::getzone(size_t row_start,size_t line_start,size_t row_end,size_t line_end) const
{
	if(row_start==0||row_end==0||line_start==0||line_start==0) throw bad_oper("getzone出错:行列标不能为0!");
	if(row_start>row_end) std::swap(row_start,row_end);
	if(line_start>line_end) std::swap(line_start,line_end);
	if (row_end>row_count) throw bad_oper("getzone出错:行标越界!");
	if (line_end>line_count) throw bad_oper("getzone出错:列标越界!");

	Cmatrix<Type> tmp(row_end-row_start+1,line_end-line_start+1);
	for(size_t i=1;i<=tmp.row();i++)
		for (size_t j=1;j<=tmp.line();j++)
			tmp(i,j)=(*this)(row_start+i-1,line_start+j-1);
	return tmp;
};

template<class Type>
inline void Cmatrix<Type>::Lmax(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg,size_t end) const
{
	if(beg==0||end==0) throw bad_oper("Lmax出错:行标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>row_count) throw bad_oper("Lmax出错:行标越界!");
	rst.resize(1,line_count);
	pos.resize(1,line_count);
	for (size_t i=1;i<=line_count;i++)
	{
		Type value=(*this)(beg,i);
		size_t _pos=beg;
		for (size_t j=beg+1;j<=end;j++)
		{
			if ((*this)(j,i)>value)
			{
				value=(*this)(j,i);
				_pos=j;
			}
		}
		rst(1,i)=value;
		pos(1,i)=_pos;
	}
};

template<class Type>
inline void Cmatrix<Type>::Lmaxa(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg,size_t end) const
{
	if(beg==0||end==0) throw bad_oper("Lmaxa出错:行标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>row_count) throw bad_oper("Lmaxa出错:行标越界!");
	rst.resize(1,line_count);
	pos.resize(1,line_count);
	for (size_t i=1;i<=line_count;i++)
	{
		Type value=(*this)(beg,i);
		size_t _pos=beg;
		for (size_t j=beg+1;j<=end;j++)
		{
			if (aufunc::abs((*this)(j,i))>aufunc::abs(value))
			{
				value=(*this)(j,i);
				_pos=j;
			}
		}
		rst(1,i)=value;
		pos(1,i)=_pos;
	}
};


template<class Type>
inline void Cmatrix<Type>::Rmax(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg,size_t end) const
{
	if(beg==0||end==0) throw bad_oper("Rmax出错:列标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>line_count) throw bad_oper("Rmax出错:列标越界!");
	rst.resize(row_count,1);
	pos.resize(row_count,1);
	for (size_t i=1;i<=row_count;i++)
	{
		Type value=(*this)(i,beg);
		size_t _pos=beg;
		for (size_t j=beg+1;j<=end;j++)
		{
			if ((*this)(i,j)>value)
			{
				value=(*this)(i,j);
				_pos=j;
			}
		}
		rst(i,1)=value;
		pos(i,1)=_pos;
	}
};

template<class Type>
inline void Cmatrix<Type>::Rmaxa(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg,size_t end) const
{
	if(beg==0||end==0) throw bad_oper("Rmaxa出错:列标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>line_count) throw bad_oper("Rmaxa出错:列标越界!");
	rst.resize(row_count,1);
	pos.resize(row_count,1);
	for (size_t i=1;i<=row_count;i++)
	{
		Type value=(*this)(i,beg);
		size_t _pos=beg;
		for (size_t j=beg+1;j<=end;j++)
		{
			if (aufunc::abs((*this)(i,j))>aufunc::abs(value))
			{
				value=(*this)(i,j);
				_pos=j;
			}
		}
		rst(i,1)=value;
		pos(i,1)=_pos;
	}
};

template<class Type>
inline void Cmatrix<Type>::Mmax(Type & rst,std::pair<size_t,size_t>& pos,Method flag) const
{
	size_t rpos=1,lpos=1;
	Type value=(*this)(1,1);
	if(flag==LINE)
	{
		for (size_t i=1;i<=line_count;i++)
		for (size_t j=1;j<=row_count;j++)
			if ((*this)(j,i)>value)
			{
				value=(*this)(j,i);
				rpos=j;
				lpos=i;
			};
	}
	else
	{
		for (size_t i=1;i<=row_count;i++)
		for (size_t j=1;j<=line_count;j++)
			if ((*this)(i,j)>value)
			{
				value=(*this)(i,j);
				rpos=i;
				lpos=j;
			};
	}
	
	pos.first=rpos;pos.second=lpos;
	rst=value;
};

template<class Type>
inline void Cmatrix<Type>::Mmaxa(Type& rst,std::pair<size_t,size_t>& pos,Method flag) const
{
	size_t rpos=1,lpos=1;
	Type value=(*this)(1,1);
	if(flag==LINE)
	{
		for (size_t i=1;i<=line_count;i++)
		for (size_t j=1;j<=row_count;j++)
			if (aufunc::abs((*this)(j,i))>aufunc::abs(value))
			{
				value=(*this)(j,i);
				rpos=j;
				lpos=i;
			};
	}
	else
	{
		for (size_t i=1;i<=row_count;i++)
		for (size_t j=1;j<=line_count;j++)
			if (aufunc::abs((*this)(i,j))>aufunc::abs(value))
			{
				value=(*this)(i,j);
				rpos=i;
				lpos=j;
			};
	}
	pos.first=rpos;pos.second=lpos;
	rst=value;
};

template<class Type>
inline void Cmatrix<Type>::Lmin(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg,size_t end) const
{
	if(beg==0||end==0) throw bad_oper("Lmin出错:行标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>row_count) throw bad_oper("Lmin出错:行标越界!");
	rst.resize(1,line_count);
	pos.resize(1,line_count);
	for (size_t i=1;i<=line_count;i++)
	{
		Type value=(*this)(beg,i);
		size_t _pos=beg;
		for (size_t j=1+beg;j<=end;j++)
		{
			if ((*this)(j,i)<value)
			{
				value=(*this)(j,i);
				_pos=j;
			}
		}
		rst(1,i)=value;
		pos(1,i)=_pos;
	}
};

template<class Type>
inline void Cmatrix<Type>::Lmina(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg,size_t end) const
{
	if(beg==0||end==0) throw bad_oper("Lmina出错:行标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>row_count) throw bad_oper("Lmina出错:行标越界!");
	rst.resize(1,line_count);
	pos.resize(1,line_count);
	for (size_t i=1;i<=line_count;i++)
	{
		Type value=(*this)(beg,i);
		size_t _pos=beg;
		for (size_t j=1+beg;j<=end;j++)
		{
			if (aufunc::abs((*this)(j,i))<aufunc::abs(value))
			{
				value=(*this)(j,i);
				_pos=j;
			}
		}
		rst(1,i)=value;
		pos(1,i)=_pos;
	}
};

template<class Type>
inline void Cmatrix<Type>::Rmin(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg,size_t end) const
{
	if(beg==0||end==0) throw bad_oper("Rmin出错:列标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>line_count) throw bad_oper("Rmin出错:列标越界!");
	rst.resize(row_count,1);
	pos.resize(row_count,1);
	for (size_t i=1;i<=row_count;i++)
	{
		Type value=(*this)(i,beg);
		size_t _pos=beg;
		for (size_t j=1+beg;j<=end;j++)
		{
			if ((*this)(i,j)<value)
			{
				value=(*this)(i,j);
				_pos=j;
			}
		}
		rst(i,1)=value;
		pos(i,1)=_pos;
	}
};

template<class Type>
inline void Cmatrix<Type>::Rmina(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg,size_t end) const
{
	if(beg==0||end==0) throw bad_oper("Rmina出错:列标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>line_count) throw bad_oper("Rmina出错:列标越界!");
	rst.resize(row_count,1);
	pos.resize(row_count,1);
	for (size_t i=1;i<=row_count;i++)
	{
		Type value=(*this)(i,beg);
		size_t _pos=beg;
		for (size_t j=1+beg;j<=end;j++)
		{
			if (aufunc::abs((*this)(i,j))<aufunc::abs(value))
			{
				value=(*this)(i,j);
				_pos=j;
			}
		}
		rst(i,1)=value;
		pos(i,1)=_pos;
	}
};

template<class Type>
inline void Cmatrix<Type>::Mmin(Type& rst,std::pair<size_t,size_t>& pos,Method flag) const
{
	size_t rpos=1,lpos=1;
	Type value=(*this)(1,1);
	if(flag==LINE)
	{
		for (size_t i=1;i<=line_count;i++)
		for (size_t j=1;j<=row_count;j++)
			if ((*this)(j,i)<value)
			{
				value=(*this)(j,i);
				rpos=j;
				lpos=i;
			};
	}
	else
	{
		for (size_t i=1;i<=row_count;i++)
		for (size_t j=1;j<=line_count;j++)
			if ((*this)(i,j)<value)
			{
				value=(*this)(i,j);
				rpos=i;
				lpos=j;
			};
	}
	
	pos.first=rpos;pos.second=lpos;
	rst=value;
};

template<class Type>
inline void Cmatrix<Type>::Mmina(Type& rst,std::pair<size_t,size_t>& pos,Method flag) const
{
	size_t rpos=1,lpos=1;
	Type value=(*this)(1,1);
	if(flag==LINE)
	{
		for (size_t i=1;i<=line_count;i++)
		for (size_t j=1;j<=row_count;j++)
			if (aufunc::abs((*this)(j,i))<aufunc::abs(value))
			{
				value=(*this)(j,i);
				rpos=j;
				lpos=i;
			};
	}
	else
	{
		for (size_t i=1;i<=row_count;i++)
		for (size_t j=1;j<=line_count;j++)
			if (aufunc::abs((*this)(i,j))<aufunc::abs(value))
			{
				value=(*this)(i,j);
				rpos=i;
				lpos=j;
			};
	}
	pos.first=rpos;pos.second=lpos;
	rst=value;
};

template<class Type>
inline std::pair<Type,size_t> Cmatrix<Type>::Lmax(const size_t _pos,size_t beg,size_t end) const
{
	if(_pos==0) throw bad_oper("Lmax出错:列标不能为0!");
	if(_pos>line_count) throw bad_oper("Lmax出错:列标越界!");
	if(beg==0||end==0) throw bad_oper("Lmax出错:行标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>row_count) throw bad_oper("Lmax出错:行标越界!");
	Type vmax=(*this)(beg,_pos);
	size_t pos=beg;
	for(size_t i=beg;i<=end;i++)
	{
		if (vmax<(*this)(i,_pos))
		{
			vmax=(*this)(i,_pos);
			pos=i;
		}
	}
	return std::make_pair(vmax,pos);
}

template<class Type>
inline std::pair<Type,size_t> Cmatrix<Type>::Lmaxa(const size_t _pos,size_t beg,size_t end) const
{
	if(_pos==0) throw bad_oper("Lmaxa出错:列标不能为0!");
	if(_pos>line_count) throw bad_oper("Lmaxa出错:列标越界!");
	if(beg==0||end==0) throw bad_oper("Lmaxa出错:行标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>row_count) throw bad_oper("Lmaxa出错:行标越界!");
	Type vmax=(*this)(beg,_pos);
	size_t pos=beg;
	for(size_t i=beg;i<=end;i++)
	{
		if (aufunc::abs(vmax)<aufunc::abs((*this)(i,_pos)))
		{
			vmax=(*this)(i,_pos);
			pos=i;
		}
	}
	return std::make_pair(vmax,pos);
}

template<class Type>
inline std::pair<Type,size_t> Cmatrix<Type>::Rmax(const size_t _pos,size_t beg,size_t end) const
{
	if(_pos==0) throw bad_oper("Rmax出错:行标不能为0!");
	if(_pos>line_count) throw bad_oper("Rmax出错:行标越界!");
	if(beg==0||end==0) throw bad_oper("Rmax出错:列标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>line_count) throw bad_oper("Rmax出错:列标越界!");
	Type vmax=(*this)(_pos,beg);
	size_t pos=beg;
	for(size_t i=beg;i<=end;i++)
	{
		if (vmax<(*this)(_pos,i))
		{
			vmax=(*this)(_pos,i);
			pos=i;
		}
	}
	return std::make_pair(vmax,pos);
}

template<class Type>
inline std::pair<Type,size_t> Cmatrix<Type>::Rmaxa(const size_t _pos,size_t beg,size_t end) const
{
	if(_pos==0) throw bad_oper("Rmaxa出错:行标不能为0!");
	if(_pos>row_count) throw bad_oper("Rmaxa出错:行标越界!");
	if(beg==0||end==0) throw bad_oper("Rmaxa出错:列标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>line_count) throw bad_oper("Rmaxa出错:列标越界!");
	Type vmax=(*this)(_pos,beg);
	size_t pos=beg;
	for(size_t i=beg;i<=end;i++)
	{
		if (aufunc::abs(vmax)<aufunc::abs((*this)(_pos,i)))
		{
			vmax=(*this)(_pos,i);
			pos=i;
		}
	}
	return std::make_pair(vmax,pos);
}

template<class Type>
inline std::pair<Type,size_t> Cmatrix<Type>::Lmin(const size_t _pos,size_t beg,size_t end) const
{
	if(_pos==0) throw bad_oper("Lmin出错:列标不能为0!");
	if(_pos>line_count) throw bad_oper("Lmin出错:列标越界!");
	if(beg==0||end==0) throw bad_oper("Lmin出错:行标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>row_count) throw bad_oper("Lmin出错:行标越界!");
	Type vmin=(*this)(beg,_pos);
	size_t pos=beg;
	for(size_t i=beg;i<=end;i++)
	{
		if (vmin>(*this)(i,_pos))
		{
			vmin=(*this)(i,_pos);
			pos=i;
		}
	}
	return std::make_pair(vmin,pos);
}

template<class Type>
inline std::pair<Type,size_t> Cmatrix<Type>::Lmina(const size_t _pos,size_t beg,size_t end) const
{
	if(_pos==0) throw bad_oper("Lmina出错:列标不能为0!");
	if(_pos>line_count) throw bad_oper("Lmina出错:列标越界!");
	if(beg==0||end==0) throw bad_oper("Lmina出错:行标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>row_count) throw bad_oper("Lmina出错:行标越界!");
	Type vmin=(*this)(beg,_pos);
	size_t pos=beg;
	for(size_t i=beg;i<=end;i++)
	{
		if (aufunc::abs(vmin)>aufunc::abs((*this)(i,_pos)))
		{
			vmin=(*this)(i,_pos);
			pos=i;
		}
	}
	return std::make_pair(vmin,pos);
}

template<class Type>
inline std::pair<Type,size_t> Cmatrix<Type>::Rmin(const size_t _pos,size_t beg,size_t end) const
{
	if(_pos==0) throw bad_oper("Rmin出错:行标不能为0!");
	if(_pos>row_count) throw bad_oper("Rmin出错:行标越界!");
	if(beg==0||end==0) throw bad_oper("Rmin出错:列标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>line_count) throw bad_oper("Rmin出错:列标越界!");
	Type vmin=(*this)(_pos,beg);
	size_t pos=beg;
	for(size_t i=beg;i<=end;i++)
	{
		if (vmin>(*this)(_pos,i))
		{
			vmin=(*this)(_pos,i);
			pos=i;
		}
	}
	return std::make_pair(vmin,pos);
}

template<class Type>
inline std::pair<Type,size_t> Cmatrix<Type>::Rmina(const size_t _pos,size_t beg,size_t end) const
{
	if(_pos==0) throw bad_oper("Rmina出错:行标不能为0!");
	if(_pos>row_count) throw bad_oper("Rmina出错:行标越界!");
	if(beg==0||end==0) throw bad_oper("Rmina出错:列标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>line_count) throw bad_oper("Rmina出错:列标越界!");
	Type vmin=(*this)(_pos,beg);
	size_t pos=beg;
	for(size_t i=beg;i<=end;i++)
	{
		if (aufunc::abs(vmin)>aufunc::abs((*this)(_pos,i)))
		{
			vmin=(*this)(_pos,i);
			pos=i;
		}
	}
	return std::make_pair(vmin,pos);
}

template<class Type>
inline void Cmatrix<Type>::Lfind(Cmatrix<size_t>&rst,Type v1,Type v2,std::pair<size_t,size_t>lpos,std::pair<size_t,size_t>rpos,Method flag) const
{
	if(lpos.first==0||lpos.second==0) throw bad_oper("Lfind出错:列标不能为0!");
	if(rpos.first==0||rpos.second==0) throw bad_oper("Lfind出错:行标不能为0!");
	if(lpos.first>lpos.second) std::swap(lpos.first,lpos.second);
	if(rpos.first>rpos.second) std::swap(rpos.first,rpos.second);
	if(lpos.second>line_count)
		throw bad_oper("Lfind出错:列标越界!");
	if(rpos.second>row_count)
		throw bad_oper("Lfind出错:行标越界!");

	if(v1>v2) std::swap(v1,v2);

	rst.resize(1,lpos.second-lpos.first+1);
	rst=0;
	
	if(flag==FW)
	{
		for(size_t i=lpos.first;i<=lpos.second;i++)
		{
			for(size_t j=rpos.first;j<=rpos.second;j++)
				if((*this)(j,i)>=v1&&(*this)(j,i)<=v2)
				{
					rst(1,i-lpos.first+1)=j;
					break;
				}
		}
	}
	else
	{
		for(size_t i=lpos.first;i<=lpos.second;i++)
		{
			for(size_t j=rpos.second;j>=rpos.first;j--)
				if((*this)(j,i)>=v1&&(*this)(j,i)<=v2)
				{
					rst(1,i-lpos.first+1)=j;
					break;
				}
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::Lfinda(Cmatrix<size_t>&rst,long double v1,long double v2,std::pair<size_t,size_t>lpos,std::pair<size_t,size_t>rpos,Method flag) const
{
	if(lpos.first==0||lpos.second==0) throw bad_oper("Lfinda出错:列标不能为0!");
	if(rpos.first==0||rpos.second==0) throw bad_oper("Lfinda出错:行标不能为0!");
	if(lpos.first>lpos.second) std::swap(lpos.first,lpos.second);
	if(rpos.first>rpos.second) std::swap(rpos.first,rpos.second);
	if(lpos.second>line_count)
		throw bad_oper("Lfinda出错:列标越界!");
	if(rpos.second>row_count)
		throw bad_oper("Lfinda出错:行标越界!");

	if(v1>v2) std::swap(v1,v2);

	rst.resize(1,lpos.second-lpos.first+1);
	rst=0;
	
	if(flag==FW)
	{
		for(size_t i=lpos.first;i<=lpos.second;i++)
		{
			for(size_t j=rpos.first;j<=rpos.second;j++)
				if(aufunc::abs((*this)(j,i))>=v1&&aufunc::abs((*this)(j,i))<=v2)
				{
					rst(1,i-lpos.first+1)=j;
					break;
				}
		}
	}
	else
	{
		for(size_t i=lpos.first;i<=lpos.second;i++)
		{
			for(size_t j=rpos.second;j>=rpos.first;j--)
				if(aufunc::abs((*this)(j,i))>=v1&&aufunc::abs((*this)(j,i))<=v2)
				{
					rst(1,i-lpos.first+1)=j;
					break;
				}
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::Rfind(Cmatrix<size_t>&rst,Type v1,Type v2,std::pair<size_t,size_t>rpos,std::pair<size_t,size_t>lpos,Method flag) const
{
	if(lpos.first==0||lpos.second==0) throw bad_oper("Rfind出错:列标不能为0!");
	if(rpos.first==0||rpos.second==0) throw bad_oper("Rfind出错:行标不能为0!");
	if(lpos.first>lpos.second) std::swap(lpos.first,lpos.second);
	if(rpos.first>rpos.second) std::swap(rpos.first,rpos.second);
	if(lpos.second>line_count)
		throw bad_oper("Rfind出错:列标越界!");
	if(rpos.second>row_count)
		throw bad_oper("Rfind出错:行标越界!");

	if(v1>v2) std::swap(v1,v2);

	rst.resize(rpos.second-rpos.first+1,1);
	rst=0;
	
	if(flag==FW)
	{
		for(size_t i=rpos.first;i<=rpos.second;i++)
		{
			for(size_t j=lpos.first;j<=lpos.second;j++)
				if((*this)(i,j)>=v1&&(*this)(i,j)<=v2)
				{
					rst(i-rpos.first+1,1)=j;
					break;
				}
		}
	}
	else
	{
		for(size_t i=rpos.first;i<=rpos.second;i++)
		{
			for(size_t j=lpos.second;j>=lpos.first;j--)
				if((*this)(i,j)>=v1&&(*this)(i,j)<=v2)
				{
					rst(i-rpos.first+1,1)=j;
					break;
				}
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::Rfinda(Cmatrix<size_t>&rst,long double v1,long double v2,std::pair<size_t,size_t>rpos,std::pair<size_t,size_t>lpos,Method flag) const
{
	if(lpos.first==0||lpos.second==0) throw bad_oper("Rfinda出错:列标不能为0!");
	if(rpos.first==0||rpos.second==0) throw bad_oper("Rfinda出错:行标不能为0!");
	if(lpos.first>lpos.second) std::swap(lpos.first,lpos.second);
	if(rpos.first>rpos.second) std::swap(rpos.first,rpos.second);
	if(lpos.second>line_count)
		throw bad_oper("Rfinda出错:列标越界!");
	if(rpos.second>row_count)
		throw bad_oper("Rfinda出错:行标越界!");

	if(v1>v2) std::swap(v1,v2);

	rst.resize(rpos.second-rpos.first+1,1);
	rst=0;
	
	if(flag==FW)
	{
		for(size_t i=rpos.first;i<=rpos.second;i++)
		{
			for(size_t j=lpos.first;j<=lpos.second;j++)
				if(aufunc::abs((*this)(i,j))>=v1&&aufunc::abs((*this)(i,j))<=v2)
				{
					rst(i-rpos.first+1,1)=j;
					break;
				}
		}
	}
	else
	{
		for(size_t i=rpos.first;i<=rpos.second;i++)
		{
			for(size_t j=lpos.second;j>=lpos.first;j--)
				if(aufunc::abs((*this)(i,j))>=v1&&aufunc::abs((*this)(i,j))<=v2)
				{
					rst(i-rpos.first+1,1)=j;
					break;
				}
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::Lfind(Cmatrix<size_t>&rst,Type v1,Method op,std::pair<size_t,size_t>lpos,std::pair<size_t,size_t>rpos,Method flag) const
{
	if(lpos.first==0||lpos.second==0) throw bad_oper("Lfind出错:列标不能为0!");
	if(rpos.first==0||rpos.second==0) throw bad_oper("Lfind出错:行标不能为0!");
	if(lpos.first>lpos.second) std::swap(lpos.first,lpos.second);
	if(rpos.first>rpos.second) std::swap(rpos.first,rpos.second);
	if(lpos.second>line_count)
		throw bad_oper("Lfind出错:列标越界!");
	if(rpos.second>row_count)
		throw bad_oper("Lfind出错:行标越界!");

	rst.resize(1,lpos.second-lpos.first+1);
	rst=0;
	
	if(flag==FW)
	{
		if(op==GT)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if((*this)(j,i)>v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==GE)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if((*this)(j,i)>=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==LT)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if((*this)(j,i)<v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==LE)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if((*this)(j,i)<=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==EQ)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if((*this)(j,i)==v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if((*this)(j,i)!=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
	}
	else
	{
		if(op==LT)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if((*this)(j,i)<v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==LE)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if((*this)(j,i)<=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==GT)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if((*this)(j,i)>v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==GE)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if((*this)(j,i)>=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==EQ)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if((*this)(j,i)==v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if((*this)(j,i)!=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::Lfinda(Cmatrix<size_t>&rst,double v1,Method op,std::pair<size_t,size_t>lpos,std::pair<size_t,size_t>rpos,Method flag=FW) const
{
	if(lpos.first==0||lpos.second==0) throw bad_oper("Lfinda出错:列标不能为0!");
	if(rpos.first==0||rpos.second==0) throw bad_oper("Lfinda出错:行标不能为0!");
	if(lpos.first>lpos.second) std::swap(lpos.first,lpos.second);
	if(rpos.first>rpos.second) std::swap(rpos.first,rpos.second);
	if(lpos.second>line_count)
		throw bad_oper("Lfinda出错:列标越界!");
	if(rpos.second>row_count)
		throw bad_oper("Lfinda出错:行标越界!");

	rst.resize(1,lpos.second-lpos.first+1);
	rst=0;
	
	if(flag==FW)
	{
		if(op==GT)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if(aufunc::abs((*this)(j,i))>v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==GE)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if(aufunc::abs((*this)(j,i))>=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==LT)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if(aufunc::abs((*this)(j,i))<v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==LE)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if(aufunc::abs((*this)(j,i))<=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==EQ)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if(aufunc::abs((*this)(j,i))==v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.first;j<=rpos.second;j++)
					if(aufunc::abs((*this)(j,i))!=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
	}
	else
	{
		if(op==LT)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if(aufunc::abs((*this)(j,i))<v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==LE)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if(aufunc::abs((*this)(j,i))<=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==GT)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if(aufunc::abs((*this)(j,i))>v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==GE)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if(aufunc::abs((*this)(j,i))>=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else if(op==EQ)
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if(aufunc::abs((*this)(j,i))==v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
		else
		{
			for(size_t i=lpos.first;i<=lpos.second;i++)
			{
				for(size_t j=rpos.second;j>=rpos.first;j--)
					if(aufunc::abs((*this)(j,i))!=v1)
					{
						rst(1,i-lpos.first+1)=j;
						break;
					}
			}
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::Rfind(Cmatrix<size_t>&rst,Type v1,Method op,std::pair<size_t,size_t>rpos,std::pair<size_t,size_t>lpos,Method flag) const
{
	if(lpos.first==0||lpos.second==0) throw bad_oper("Rfind出错:列标不能为0!");
	if(rpos.first==0||rpos.second==0) throw bad_oper("Rfind出错:行标不能为0!");
	if(lpos.first>lpos.second) std::swap(lpos.first,lpos.second);
	if(rpos.first>rpos.second) std::swap(rpos.first,rpos.second);
	if(lpos.second>line_count)
		throw bad_oper("Rfind出错:列标越界!");
	if(rpos.second>row_count)
		throw bad_oper("Rfind出错:行标越界!");

	rst.resize(rpos.second-rpos.first+1,1);
	rst=0;
	
	if(flag==FW)
	{
		if(op==GT)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if((*this)(i,j)>v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==GE)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if((*this)(i,j)>=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==LT)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if((*this)(i,j)<v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==LE)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if((*this)(i,j)<=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==EQ)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if((*this)(i,j)==v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if((*this)(i,j)!=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
	}
	else
	{
		if(op==LT)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if((*this)(i,j)<v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==LE)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if((*this)(i,j)<=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==GT)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if((*this)(i,j)>v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==GE)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if((*this)(i,j)>=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==EQ)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if((*this)(i,j)==v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if((*this)(i,j)!=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::Rfinda(Cmatrix<size_t>&rst,long double v1,Method op,std::pair<size_t,size_t>rpos,std::pair<size_t,size_t>lpos,Method flag) const
{
	if(lpos.first==0||lpos.second==0) throw bad_oper("Rfinda出错:列标不能为0!");
	if(rpos.first==0||rpos.second==0) throw bad_oper("Rfinda出错:行标不能为0!");
	if(lpos.first>lpos.second) std::swap(lpos.first,lpos.second);
	if(rpos.first>rpos.second) std::swap(rpos.first,rpos.second);
	if(lpos.second>line_count)
		throw bad_oper("Rfinda出错:列标越界!");
	if(rpos.second>row_count)
		throw bad_oper("Rfinda出错:行标越界!");

	rst.resize(rpos.second-rpos.first+1,1);
	rst=0;
	
	if(flag==FW)
	{
		if(op==GT)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if(aufunc::abs((*this)(i,j))>v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==GE)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if(aufunc::abs((*this)(i,j))>=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==LT)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if(aufunc::abs((*this)(i,j))<v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==LE)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if(aufunc::abs((*this)(i,j))<=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==EQ)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if(aufunc::abs((*this)(i,j))==v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.first;j<=lpos.second;j++)
					if(aufunc::abs((*this)(i,j))!=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
	}
	else
	{
		if(op==LT)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if(aufunc::abs((*this)(i,j))<v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==LE)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if(aufunc::abs((*this)(i,j))<=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==GT)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if(aufunc::abs((*this)(i,j))>v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==GE)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if(aufunc::abs((*this)(i,j))>=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else if(op==EQ)
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if(aufunc::abs((*this)(i,j))==v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
		else
		{
			for(size_t i=rpos.first;i<=rpos.second;i++)
			{
				for(size_t j=lpos.second;j>=lpos.first;j--)
					if(aufunc::abs((*this)(i,j))!=v1)
					{
						rst(i-rpos.first+1,1)=j;
						break;
					}
			}
		}
	}
}

template<class Type>
inline std::pair<size_t,size_t> Cmatrix<Type>::Mfind(Type v1,Type v2,Method flag,Method flag2) const
{
	if(v1>v2) std::swap(v1,v2);
	if(flag2==FW)
	{
		if(flag==LINE)
		{
			for(size_t i=1;i<=line_count;++i)
				for(size_t j=1;j<=row_count;++j)
					if((*this)(j,i)>=v1&&(*this)(j,i)<=v2) return std::make_pair(j,i);
		}
		else
		{
			for(size_t i=1;i<=row_count;++i)
				for(size_t j=1;j<=line_count;++j)
					if((*this)(i,j)>=v1&&(*this)(i,j)<=v2) return std::make_pair(i,j);
		}
	}
	else
	{
		if(flag==LINE)
		{
			for(size_t i=line_count;i>=1;--i)
				for(size_t j=row_count;j>=1;--j)
					if((*this)(j,i)>=v1&&(*this)(j,i)<=v2) return std::make_pair(j,i);
		}
		else
		{
			for(size_t i=row_count;i>=1;--i)
				for(size_t j=line_count;j>=1;--j)
					if((*this)(i,j)>=v1&&(*this)(i,j)<=v2) return std::make_pair(i,j);
		}
	}
	return std::make_pair(0,0);
}

template<class Type>
inline std::pair<size_t,size_t> Cmatrix<Type>::Mfind(Type value,Method flag1,Method flag2,Method flag3) const
{
	if(flag3==FW)
	{
		if(flag2==LINE)
		{
			if(flag1==GT)
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if((*this)(j,i)>value) return std::make_pair(j,i);
			}
			else if(flag1==LT)
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if((*this)(j,i)<value) return std::make_pair(j,i);
			}
			else if(flag1==GE)
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if((*this)(j,i)>=value) return std::make_pair(j,i);
			}
			else if(flag1==LE)
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if((*this)(j,i)<=value) return std::make_pair(j,i);
			}
			else if(flag1==EQ)
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if((*this)(j,i)==value) return std::make_pair(j,i);
			}
			else
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if((*this)(j,i)!=value) return std::make_pair(j,i);
			}
		}
		else
		{
			if(flag1==GT)
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if((*this)(i,j)>value) return std::make_pair(i,j);
			}
			else if(flag1==LT)
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if((*this)(i,j)<value) return std::make_pair(i,j);
			}
			else if(flag1==GE)
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if((*this)(i,j)>=value) return std::make_pair(i,j);
			}
			else if(flag1==LE)
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if((*this)(i,j)<=value) return std::make_pair(i,j);;
			}
			else if(flag1==EQ)
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if((*this)(i,j)==value) return std::make_pair(i,j);
			}
			else
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if((*this)(i,j)!=value) return std::make_pair(i,j);
			}
		}
	}
	else
	{
		if(flag2==LINE)
		{
			if(flag1==GT)
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=row_count;j>=1;--j)
						if((*this)(j,i)>value) return std::make_pair(j,i);
			}
			else if(flag1==LT)
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=row_count;j>=1;--j)
						if((*this)(j,i)<value) return std::make_pair(j,i);
			}
			else if(flag1==GE)
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=row_count;j>=1;--j)
						if((*this)(j,i)>=value) return std::make_pair(j,i);
			}
			else if(flag1==LE)
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if((*this)(j,i)<=value) return std::make_pair(j,i);
			}
			else if(flag1==EQ)
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=row_count;j>=1;--j)
						if((*this)(j,i)==value) return std::make_pair(j,i);
			}
			else
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=row_count;j>=1;--j)
						if((*this)(j,i)!=value) return std::make_pair(j,i);
			}
		}
		else
		{
			if(flag1==GT)
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if((*this)(i,j)>value) return std::make_pair(i,j);
			}
			else if(flag1==LT)
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if((*this)(i,j)<value) return std::make_pair(i,j);
			}
			else if(flag1==GE)
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if((*this)(i,j)>=value) return std::make_pair(i,j);
			}
			else if(flag1==LE)
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if((*this)(i,j)<=value) return std::make_pair(i,j);;
			}
			else if(flag1==EQ)
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if((*this)(i,j)==value) return std::make_pair(i,j);
			}
			else
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if((*this)(i,j)!=value) return std::make_pair(i,j);
			}
		}
	}
	return std::make_pair(0,0);
}

template<class Type>
inline std::pair<size_t,size_t> Cmatrix<Type>::Mfinda(long double v1,long double v2,Method flag,Method flag2) const
{
	if(v1>v2) std::swap(v1,v2);
	if(flag2==FW)
	{
		if(flag==LINE)
		{
			for(size_t i=1;i<=line_count;++i)
				for(size_t j=1;j<=row_count;++j)
					if(aufunc::abs((*this)(j,i))>=v1&&aufunc::abs((*this)(j,i))<=v2) return std::make_pair(j,i);
		}
		else
		{
			for(size_t i=1;i<=row_count;++i)
				for(size_t j=1;j<=line_count;++j)
					if(aufunc::abs((*this)(i,j))>=v1&&aufunc::abs((*this)(i,j))<=v2) return std::make_pair(i,j);
		}
	}
	else
	{
		if(flag==LINE)
		{
			for(size_t i=line_count;i>=1;--i)
				for(size_t j=row_count;j>=1;--j)
					if(aufunc::abs((*this)(i,j))>=v1&&aufunc::abs((*this)(i,j))<=v2) return std::make_pair(j,i);
		}
		else
		{
			for(size_t i=row_count;i>=1;--i)
				for(size_t j=line_count;j>=1;--j)
					if(aufunc::abs((*this)(i,j))>=v1&&aufunc::abs((*this)(i,j))<=v2) return std::make_pair(i,j);
		}
	}
	return std::make_pair(0,0);
}

template<class Type>
inline std::pair<size_t,size_t> Cmatrix<Type>::Mfinda(long double value,Method flag1,Method flag2,Method flag3) const
{
	if(flag3==FW)
	{
		if(flag2==LINE)
		{
			if(flag1==GT)
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if(aufunc::abs((*this)(i,j))>value) return std::make_pair(j,i);
			}
			else if(flag1==LT)
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if(aufunc::abs((*this)(i,j))<value) return std::make_pair(j,i);
			}
			else if(flag1==GE)
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if(aufunc::abs((*this)(i,j))>=value) return std::make_pair(j,i);
			}
			else if(flag1==LE)
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if(aufunc::abs((*this)(i,j))<=value) return std::make_pair(j,i);
			}
			else if(flag1==EQ)
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if(aufunc::abs((*this)(i,j))==value) return std::make_pair(j,i);
			}
			else
			{
				for(size_t i=1;i<=line_count;++i)
					for(size_t j=1;j<=row_count;++j)
						if(aufunc::abs((*this)(i,j))!=value) return std::make_pair(j,i);
			}
		}
		else
		{
			if(flag1==GT)
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if(aufunc::abs((*this)(i,j))>value) return std::make_pair(i,j);
			}
			else if(flag1==LT)
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if(aufunc::abs((*this)(i,j))<value) return std::make_pair(i,j);
			}
			else if(flag1==GE)
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if(aufunc::abs((*this)(i,j))>=value) return std::make_pair(i,j);
			}
			else if(flag1==LE)
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if(aufunc::abs((*this)(i,j))<=value) return std::make_pair(i,j);;
			}
			else if(flag1==EQ)
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if(aufunc::abs((*this)(i,j))==value) return std::make_pair(i,j);
			}
			else
			{
				for(size_t i=1;i<=row_count;++i)
					for(size_t j=1;j<=line_count;++j)
						if(aufunc::abs((*this)(i,j))!=value) return std::make_pair(i,j);
			}
		}
	}
	else
	{
		if(flag2==LINE)
		{
			if(flag1==GT)
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=row_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))>value) return std::make_pair(j,i);
			}
			else if(flag1==LT)
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=row_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))<value) return std::make_pair(j,i);
			}
			else if(flag1==GE)
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=row_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))>=value) return std::make_pair(j,i);
			}
			else if(flag1==LE)
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))<=value) return std::make_pair(j,i);
			}
			else if(flag1==EQ)
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=row_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))==value) return std::make_pair(j,i);
			}
			else
			{
				for(size_t i=line_count;i>=1;--i)
					for(size_t j=row_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))!=value) return std::make_pair(j,i);
			}
		}
		else
		{
			if(flag1==GT)
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))>value) return std::make_pair(i,j);
			}
			else if(flag1==LT)
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))<value) return std::make_pair(i,j);
			}
			else if(flag1==GE)
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))>=value) return std::make_pair(i,j);
			}
			else if(flag1==LE)
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))<=value) return std::make_pair(i,j);;
			}
			else if(flag1==EQ)
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))==value) return std::make_pair(i,j);
			}
			else
			{
				for(size_t i=row_count;i>=1;--i)
					for(size_t j=line_count;j>=1;--j)
						if(aufunc::abs((*this)(i,j))!=value) return std::make_pair(i,j);
			}
		}
	}
	return std::make_pair(0,0);
}

template<class Type>
template<class InputIterator,class U>
inline Type Cmatrix<Type>::interp(InputIterator lit,InputIterator rit,const std::pair<U,U> cord) const
{
	size_t row1=0,line1=0,row2=0,line2=0;
	Type rtn;
	InputIterator tmp=lit;
	bool flag;
	if(line_count==1)
	{
		if(cord.first!=*lit) throw bad_oper("interp出错:列坐标超出范围!");
		else line1=line2=1;
	}
	else
	{
		++tmp;
		if(*lit==*tmp) throw bad_oper("interp出错:列坐标向量非单调!");
		else if(*lit<*tmp) flag=true;
		else flag=false;
		if(flag)
		{
			if(cord.first<*lit) throw bad_oper("interp出错:插值列坐标超出范围!");
			for (line2=2;line2<=line_count;line2++)
			{
				lit++;
				if(*lit>=cord.first) break;
			}
		}
		else
		{
			if(cord.first>*lit) throw bad_oper("interp出错:插值列坐标超出范围!");
			for (line2=2;line2<=line_count;line2++)
			{
				lit++;
				if(*lit<=cord.first) break;		
			}
		}
		if(line2>line_count) throw bad_oper("interp出错:插值列坐标超出范围或原坐标非单调!");
		else line1=line2-1;
	}
	
	if(row_count==1)
	{
		if(cord.second!=*rit) throw bad_oper("interp出错:行坐标超出范围!");
		else row1=row2=1;
	}
	else
	{
		tmp=rit;
		++tmp;
		if(*rit==*tmp) throw bad_oper("interp出错:行坐标向量非单调!");
		else if(*rit<*tmp) flag=true;
		else flag=false;
		if(flag)
		{
			if(cord.second<*rit) throw bad_oper("interp出错:插值行坐标超出范围!");
			for (row2=2;row2<=row_count;row2++)
			{
				rit++;
				if(*rit>=cord.second) break;
			}
		}
		else
		{
			if(cord.second>*rit) throw bad_oper("interp出错:插值行坐标超出范围!");
			for (row2=2;row2<=row_count;row2++)
			{
				rit++;
				if(*rit<=cord.second) break;		
			}
		}
		if(row2>row_count) throw bad_oper("interp出错:插值行坐标超出范围或原坐标非单调!");
		else row1=row2-1;
	}

	const Type lt=(*this)(row1,line1),rt=(*this)(row1,line2),lb=(*this)(row2,line1),rb=(*this)(row2,line2);
	Type ct,cb;
	if(line1==line2) {ct=lt;cb=lb;}
	else
	{
		tmp=lit;lit--;
		ct=lt+(rt-lt)/Type(*tmp-*lit)*Type(cord.first-*lit);
		cb=lb+(rb-lb)/Type(*tmp-*lit)*Type(cord.first-*lit);
	}

	if(row2==row1) rtn=ct;
	else
	{
		tmp=rit;rit--;
		rtn=ct+(cb-ct)/Type(*tmp-*rit)*Type(cord.second-*rit);
	}
	return rtn;
}

template<class Type>
template<class InputIterator,class U>
inline void Cmatrix<Type>::interp(InputIterator lit,InputIterator rit,const Cmatrix<std::pair<U,U> >&cord,Cmatrix<Type>& rtn) const
{
	rtn.resize(cord.row(),cord.line());
	for(size_t i=1;i<=cord.row();i++)
		for(size_t j=1;j<=cord.line();j++)
			rtn(i,j)=interp(lit,rit,cord(i,j));
}

template<class Type>
template<class U>
inline Type Cmatrix<Type>::interp(const Cmatrix<std::pair<U,U> > &cord,const std::pair<U,U>pos) const
{
	if(row_count!=cord.row()||line_count!=cord.line()) throw bad_oper("interp出错:坐标行列数不对应!");
	if(row_count<=1||line_count<=1) throw bad_oper("interp出错:矩阵至少应为2*2阶!");
	size_t i,j;
	for(i=1;i<row_count;i++)
	{
		std::pair<U,U> quad[4];
		for(j=1;j<line_count;j++)
		{
			quad[0]=cord(i,j);
			quad[1]=cord(i+1,j);
			quad[2]=cord(i+1,j+1);
			quad[3]=cord(i,j+1);
			if(aufunc::bInArea(quad,pos)) break;
		}
		if(j<line_count) break;
	}
	if(i==row_count) throw bad_oper("interp出错:坐标超出范围!");
	std::pair<U,U> xy[4]={cord(i,j),cord(i,j+1),cord(i+1,j+1),cord(i+1,j)},tmp[3];
	Type z[4]={(*this)(i,j),(*this)(i,j+1),(*this)(i+1,j+1),(*this)(i+1,j)};
	Cmatrix<Type> _z(3,1);
	Cmatrix<U> _xy(3,3);
	Type rtn=0;
	size_t num=0;
	for(i=0;i<4;i++)
	{
		size_t k=0;
		for(j=0;j<4;j++)
		{
			if(j==i) continue;
			tmp[k]=xy[j];
			k++;
		}
		if(!aufunc::bInAreaT(tmp,pos)) continue;
		k=0;
		num++;

		for(j=0;j<4;j++)
		{
			if(j==i) continue;
			_z(k+1,1)=z[j];
			_xy(k+1,1)=xy[j].first;_xy(k+1,2)=xy[j].second;_xy(k+1,3)=1;
			k++;
		}
		_xy.invered();
		_z.Rmutiply(_xy);
		rtn=rtn+_z(1,1)*pos.first+_z(2,1)*pos.second+_z(3,1);
	}
	rtn/=num;
	return rtn;
}

template<class Type>
template<class U>
inline void Cmatrix<Type>::interp(const Cmatrix<std::pair<U,U> > &cord,const Cmatrix<std::pair<U,U> >pos,Cmatrix<Type>&rst) const
{
	rst.resize(pos.row(),pos.line());
	for(size_t i=1;i<=pos.row();i++)
		for(size_t j=1;j<=pos.line();j++)
			rtn(i,j)=interp(cord,pos(i,j));
}

template<class Type>
template<class U>
inline void Cmatrix<Type>::real(Cmatrix<U>&rst) const
{
	rst.resize(row_count,line_count);
	for(size_t i=1;i<=row_count;i++)
		for(size_t j=1;j<=line_count;j++)
			rst(i,j)=aufunc::real((*this)(i,j));
}

template<class Type>
template<class U>
inline void Cmatrix<Type>::imag(Cmatrix<U>&rst) const
{
	rst.resize(row_count,line_count);
	for(size_t i=1;i<=row_count;i++)
		for(size_t j=1;j<=line_count;j++)
			rst(i,j)=aufunc::imag((*this)(i,j));
}

template<class Type>
inline Cmatrix<long double> Cmatrix<Type>::real() const
{
	Cmatrix<long double> tmp;
	real(tmp);
	return tmp;
}

template<class Type>
inline Cmatrix<long double> Cmatrix<Type>::imag() const
{
	Cmatrix<long double> tmp;
	imag(tmp);
	return tmp;
}

//变动型矩阵变换************************************************
template<class Type>
inline Type & Cmatrix<Type>::operator ()(size_type row,size_type line)
{
	if(row==0)throw bad_oper("operator()出错:行标不能为0!");
	if(line==0)throw bad_oper("operator()出错:列标不能为0!");
	return matrix[line-1][row-1];
};

template<class Type>
inline void Cmatrix<Type>::resize(const size_type row,const size_type line,Type elem=Type())
{
	matrix.resize(line);
	for(size_t i=line_count;i<line;++i)
		matrix[i].resize(row,elem);
	if(row==row_count) 
	{
		line_count=line;
		return;
	}
	for(size_t i=0;i<std::min(line,line_count);++i)
		matrix[i].resize(row,elem);
	row_count=row;
	line_count=line;
};

template<class Type>
inline void Cmatrix<Type>::Lresize(size_type line,Type elem=Type())
{
	if (line==line_count) return ;
	matrix.resize(line);
	for (size_type i=line_count;i<line;matrix[i++].resize(row_count,elem));
	line_count=line;
};

template<class Type>
inline void Cmatrix<Type>::Rresize(size_type row,Type elem=Type())
{
	if (row==row_count) return ;
	for (size_type i=0;i<line_count;matrix[i++].resize(row,elem));
	row_count=row;
};

template <class Type>
template <class	U>
inline void Cmatrix<Type>::combine(const Cmatrix<U> &orig,Method flag=LINE)
{
	const size_t row=row_count;
	const size_t line=line_count;
	if(flag==LINE)
	{
		if (row_count!=orig.row())
		throw bad_oper("combine出错:两矩阵行数不等!");
		Lresize(line_count+orig.line());
		for (size_t i=1;i<=row_count;i++)
			for (size_t j=line+1;j<=line_count;j++)
				(*this)(i,j)=orig(i,j-line);
	}
	else
	{
		if (line_count!=orig.line())
		throw bad_oper("combine出错:两矩阵列数不等!");
		Rresize(row_count+orig.row());
		for (size_t i=1+row;i<=row_count;i++)
			for (size_t j=1;j<=line_count;j++)
				(*this)(i,j)=orig(i-row,j);
	}
};

template<class Type>
inline void Cmatrix<Type>::row(size_t start,size_t end)
{
	if(start==0||end==0) throw bad_oper("row出错:行标不能为0");
	if(start>end) std::swap(start,end);
	if (end>row_count) throw bad_oper("row出错:行标越界!");
	for(size_t i=1;i<=end-start+1;i++)
		for (size_t j=1;j<=line_count;j++)
			(*this)(i,j)=(*this)(start+i-1,j);
	Rresize(end-start+1);
}

template<class Type>
inline void Cmatrix<Type>::line(size_t start,size_t end)
{
	if(start==0||end==0) throw bad_oper("line出错:列标不能为0");
	if(start>end) std::swap(start,end);
	if (end>line_count) throw bad_oper("line出错:列标越界!");
	for(size_t i=1;i<=row_count;i++)
		for (size_t j=1;j<=end-start+1;j++)
			(*this)(i,j)=(*this)(i,start+j-1);
	Lresize(end-start+1);
}

template<class Type>
inline void Cmatrix<Type>::zone(size_t row_start,size_t line_start,size_t row_end,size_t line_end)
{
	if(row_start==0||row_end==0||line_start==0||line_end==0) throw bad_oper("zone出错:行列标不能为0!");
	if(row_start>row_end) std::swap(row_start,row_end);
	if(line_start>line_end) std::swap(line_start,line_end);
	if (row_end>row_count) throw bad_oper("zone出错:行标越界!");
	if (line_end>line_count) throw bad_oper("zone出错:列标越界!");

	for(size_t i=1;i<=row_end-row_start+1;i++)
		for (size_t j=1;j<=line_end-line_start+1;j++)
			(*this)(i,j)=(*this)(row_start+i-1,line_start+j-1);
	resize(row_end-row_start+1,line_end-line_start+1);
};


template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Linsert(const size_t pos,const size_t num,InputIterator beg,InputIterator end,const Type&elem)
{
	if(pos==0) throw bad_oper("Linsert出错:列标不能为0!");
	if (pos>line_count+1) throw bad_oper("Linsert出错:列标越界!");
	std::vector<Type> tmp(row_count,elem);
	std::vector<vector<Type> >::iterator it;
	it=matrix.begin()+(pos-1);
	matrix.insert(it,num,tmp);
	line_count+=num;
	for (size_t i=pos;i<pos+num;i++)
	{
		if(beg==end) break;
		for (size_t j=1;j<=row_count;j++)
		{
			if (beg!=end) (*this)(j,i)=*beg++;
			else break;
		}
	}
};

template<class Type>
template<class U>
inline void Cmatrix<Type>::Linsert(const size_t pos,const size_t num,const U* ptr,size_t num2,const Type&elem)
{
	if(pos==0) throw bad_oper("Linsert出错:列标不能为0!");
	if (pos>line_count+1) throw bad_oper("Linsert出错:列标越界!");
	std::vector<Type> tmp(row_count);
	std::vector<vector<Type> >::iterator it;
	it=matrix.begin()+(pos-1);
	matrix.insert(it,num,tmp);
	line_count+=num;
	for (size_t i=pos;i<pos+num;i++)
	{
		if(num2==0) break;
		for (size_t j=1;j<=row_count;j++)
		{
			if (num2!=0) (*this)(j,i)=*ptr++;
			else break;
			num2--;
		}
	}
};

template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Linsert(const size_t pos,const size_t num,InputIterator beg)
{
	if(pos==0) throw bad_oper("Linsert出错:列标不能为0!");
	if (pos>line_count+1) throw bad_oper("Linsert出错:列标越界!");
	std::vector<Type> tmp(row_count);
	std::vector<vector<Type> >::iterator it;
	it=matrix.begin()+(pos-1);
	matrix.insert(it,num,tmp);
	line_count+=num;
	for (size_t i=pos;i<pos+num;i++)
		for (size_t j=1;j<=row_count;j++)
			(*this)(j,i)=*beg++;
};


template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Rinsert(const size_t pos,const size_t num,InputIterator beg,InputIterator end,const Type&elem)
{
	if(pos==0) throw bad_oper("Rinsert出错:行标不能为0!");
	if (pos>row_count+1) throw bad_oper("Rinsert出错:行标越界!");
	for (size_t i=1;i<=num;i++)
		for (size_t j=1;j<=line_count;j++)
			if (beg==end) {matrix[j-1].insert(matrix[j-1].begin()+(pos+i-2),Type());continue;}
			else matrix[j-1].insert(matrix[j-1].begin()+(pos+i-2),*beg++);
	row_count+=num;
}

template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Rinsert(const size_t pos,const size_t num,InputIterator beg)
{
	if(pos==0) throw bad_oper("Rinsert出错:行标不能为0!");
	if (pos>row_count+1) throw bad_oper("Rinsert出错:行标越界!");
	for (size_t i=1;i<=num;i++)
		for (size_t j=1;j<=line_count;j++)
			matrix[j-1].insert(matrix[j-1].begin()+(pos+i-2),*beg++);
	row_count+=num;
}

template<class Type>
template<class U>
inline void Cmatrix<Type>::Rinsert(const size_t pos,const size_t num,const U* ptr,size_t num2,const Type&elem)
{
	if(pos==0) throw bad_oper("Rinsert出错:行标不能为0!");
	if (pos>row_count+1) throw bad_oper("Rinsert出错:行标越界!");
	for (size_t i=1;i<=num;i++)
		for (size_t j=1;j<=line_count;j++)
		{
			if (num2==0) {matrix[j-1].insert(matrix[j-1].begin()+(pos+i-2),Type());continue;}
			else matrix[j-1].insert(matrix[j-1].begin()+(pos+i-2),*ptr++);
			num2--;
		}
	row_count+=num;
};

template<class Type>
inline void Cmatrix<Type>::Lerase(const size_t pos)
{
	if(pos==0) throw bad_oper("Lerase出错:列标不能为0!");
	if(pos>line_count) throw bad_oper("Lerase出错:列标越界!")
	matrix.erase(matrix.begin()+(pos-1));
	--line_count;
};

template<class Type>
inline void Cmatrix<Type>::Lerase(size_t pos1,size_t pos2)
{
	if(pos1==0||pos2==0) throw bad_oper("Lerase出错:列标不能为0!");
	if(pos1>pos2) std::swap(pos1,pos2);
	if(pos2>line_count) throw bad_oper("Lerase出错:列标越界!");
	matrix.erase(matrix.begin()+(pos1-1),matrix.begin()+pos2);
	line_count-=(pos2-pos1+1);
};

template<class Type>
inline void Cmatrix<Type>::Rerase(const size_t pos)
{
	if(pos==0) throw bad_oper("Rerase出错:行标不能为0!");
	if(pos>row_count) throw bad_oper("Rerase出错:矩阵行标越界!");
	--row_count;
	for (size_t i=1;i<=line_count;i++)
		matrix[i-1].erase(matrix[i-1].begin()+(pos-1));
};

template<class Type>
inline void Cmatrix<Type>::Rerase(size_t pos1,size_t pos2)
{
	if(pos1==0||pos2==0) throw bad_oper("Rerase出错:列标不能为0!");
	if(pos1>pos2) std::swap(pos1,pos2);
	if(pos2>row_count) throw bad_oper("Rerase出错:列标越界!");
	for (size_t i=1;i<=line_count;i++)
		matrix[i-1].erase(matrix[i-1].begin()+(pos1-1),matrix[i-1].begin()+pos2);
	row_count-=(pos2-pos1+1);
};

template<class Type>
inline void Cmatrix<Type>::Radd(const size_t pos1,const size_t pos2,const size_t pos_r,const Type factor1,const Type factor2)
{
	if(pos1==0||pos2==0||pos_r==0) throw bad_oper("Radd出错:行标不能为0!");
	if(pos1>row_count||pos2>row_count||pos_r>row_count) throw bad_oper("Radd出错:行标越界!");
	for (size_t i=1;i<=line_count;i++)
		(*this)(pos_r,i)=(*this)(pos1,i)*factor1+(*this)(pos2,i)*factor2;
}

template<class Type>
inline void Cmatrix<Type>::Ladd(const size_t pos1,const size_t pos2,const size_t pos_r,const Type factor1=Type(1),const Type factor2=Type(1))
{
	if(pos1==0||pos2==0||pos_r==0) throw bad_oper("Ladd出错:列标不能为0!");
	if(pos1>line_count||pos2>line_count||pos_r>line_count) throw bad_oper("Ladd出错:列标越界!");
	for (size_t i=1;i<=row_count;i++)
		(*this)(i,pos_r)=(*this)(i,pos1)*factor1+(*this)(i,pos2)*factor2;
}

template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Lpush_back(InputIterator beg)
{
	InputIterator end(beg);
	std::advance(end,row_count);
	vector<Type> tmp(beg,end);
	matrix.push_back(tmp);
	++line_count;
}

template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Lpush_back(InputIterator beg,InputIterator end,const Type & elem)
{
	vector<Type> tmp(beg,end);
	tmp.resize(row_count,elem);
	matrix.push_back(tmp);
	++line_count;
}

template<class Type>
template<class U>
inline void Cmatrix<Type>::Lpush_back(const U* beg,size_t num,const Type& elem)
{
	std::vector<Type> tmp(beg,beg+num);
	tmp.resize(row_count,elem);
	matrix.push_back(tmp);
	++line_count;

	
}

template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Rpush_back(InputIterator beg)
{
	for(size_t i=0;i<line_count;++i)
			matrix[i].push_back(*beg++);
	++row_count;
}

template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Rpush_back(InputIterator beg,InputIterator end,const Type & elem)
{
	for(size_t i=0;i<line_count;++i)
			if(beg!=end) matrix[i].push_back(*beg++);
			else  matrix[i].push_back(elem);
	++row_count;
}

template<class Type>
template<class U>
inline void Cmatrix<Type>::Rpush_back(const U* beg,size_t num,const Type& elem)
{
	for(size_t i=0;i<line_count;++i)
	{
		if(num!=0) {matrix[i].push_back(*beg++);num--;}
		else  
			matrix[i].push_back(elem);
	}
	++row_count;
}

template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Lpush_up(InputIterator beg)
{
	InputIterator end(beg);
	std::advance(end,row_count);
	vector<Type> tmp(beg,end);
	matrix.insert(matrix.begin(),tmp);
	++line_count;
}

template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Lpush_up(InputIterator beg,InputIterator end,const Type& elem)
{
	vector<Type> tmp(beg,end);
	tmp.resize(row_count,elem);
	matrix.insert(matrix.begin(),tmp);
	++line_count;
}

template<class Type>
template<class U>
inline void Cmatrix<Type>::Lpush_up(const U* beg,size_t num,const Type& elem)
{
	vector<Type> tmp(beg,beg+num);
	tmp.resize(row_count,elem);
	matrix.insert(matrix.begin(),tmp);
	++line_count;
}

template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Rpush_up(InputIterator beg)
{
	for(size_t i=0;i<line_count;++i)
		matrix[i].insert(matrix[i].begin(),*beg++);
	++row_count;
}

template<class Type>
template<class InputIterator>
inline void Cmatrix<Type>::Rpush_up(InputIterator beg,InputIterator end,const Type& elem)
{
	for(size_t i=0;i<line_count;++i)
		if(beg!=end) matrix[i].insert(matrix[i].begin(),*beg++);
		else  matrix[i].insert(matrix[i].begin(),elem);
	++row_count;
}

template<class Type>
template<class U>
inline void Cmatrix<Type>::Rpush_up(const U* beg,size_t num,const Type& elem)
{
	for(size_t i=0;i<line_count;++i)
		if(num!=0) {matrix[i].insert(matrix[i].begin(),*beg++);num--;}
		else  matrix[i].insert(matrix[i].begin(),elem);
	++row_count;
}

template<class Type>
inline void Cmatrix<Type>::pop_back(Method flag=LINE)
{
	if(flag==LINE)
	{
		if(line_count==0) throw bad_oper("pop_back出错:矩阵列数为0!");
		matrix.pop_back();
		--line_count;
	}
	else
	{
		if(row_count==0) throw bad_oper("pop_back出错:矩阵行数为0!");
		for(size_t i=0;i<line_count;++i)
			matrix[i].pop_back();
		--row_count;
	}
}

template<class Type>
inline void Cmatrix<Type>::pop_up(Method flag=LINE)
{
	if(flag==LINE)
	{
		if(line_count==0) throw bad_oper("pop_back出错:矩阵列数为0!");
		matrix.erase(matrix.begin());
		--line_count;
	}
	else
	{
		if(row_count==0) throw bad_oper("pop_back出错:矩阵行数为0!");
		for(size_t i=0;i<line_count;++i)
			matrix[i].erase(matrix[i].begin());
		--row_count;
	}
}

template<class Type>
inline void Cmatrix<Type>::Mswap(Cmatrix<Type> &orig)
{
	if(this==&orig) return;
	if(line_count==0)
	{
		const size_t row=row_count;
		*this=orig;
		orig.resize(row,0);
		return;
	}
	if(row_count==0)
	{
		const size_t line=line_count;
		*this=orig;
		orig.resize(0,line);
		return;
	}
	if(orig.IsEmpty())
	{
		orig.Mswap(*this);
		return;
	}
	const size_t row1=row_count,row2=orig.row(),line1=line_count,line2=orig.line();
	const size_t row_max=std::max(row1,row2),line_max=std::max(line1,line2);
	const size_t row_min=std::min(row1,row2),line_min=std::min(line1,line2);
	for (size_t i=1;i<=row_min;i++)
		for (size_t j=1;j<=line_min;j++)
			std::swap((*this)(i,j),orig(i,j));
	if(row1<row2)
	{
		if(line1<line2)
		{
			resize(row2,line2);
		    for (size_t i=1;i<=line1;i++)
				for (size_t j=row1+1;j<=row2;j++)
					(*this)(j,i)=orig(j,i);	
			for (size_t i=line1+1;i<=line2;i++)
				for (size_t j=1;j<=row2;j++)
					(*this)(j,i)=orig(j,i);
			orig.resize(row1,line1);
		}
		else
		{
			for(size_t i=1;i<=line2;++i)
			{
				matrix[i-1].resize(row2);
				row_count=row2;
				for(size_t j=row1+1;j<=row2;j++)
				(*this)(j,i)=orig(j,i);
			}
			orig.resize(row1,line1);
		    for(size_t i=line2+1;i<=line1;++i)
				for (size_t j=1;j<=row1;j++)
					orig(j,i)=(*this)(j,i);
			Lresize(line2);
		}
	}
	else
	{
		if(line1>line2)
		{
			orig.resize(row1,line1);
		    for (size_t i=1;i<=line2;i++)
				for (size_t j=row2+1;j<=row1;j++)
					(orig)(j,i)=(*this)(j,i);	
			for (size_t i=line2+1;i<=line1;i++)
				for (size_t j=1;j<=row1;j++)
					(orig)(j,i)=(*this)(j,i);
			resize(row2,line2);
		}
		else
		{
			for(size_t i=1;i<=line1;++i)
			{
				orig.matrix[i-1].resize(row1);
				orig.row_count=row1;
				for(size_t j=row2+1;j<=row1;j++)
				orig(j,i)=(*this)(j,i);
			}
			resize(row2,line2);
		    for(size_t i=line1+1;i<=line2;++i)
				for (size_t j=1;j<=row2;j++)
					(*this)(j,i)=orig(j,i);
			orig.Lresize(line1);
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::Rswap(const size_t pos1,const size_t pos2)
{
	if(pos1==0||pos2==0) throw bad_oper("Rswap出错:行标不能为0!");
	if(pos1>row_count||pos2>row_count) throw bad_oper("Rswap出错:行标越界!");
	if(pos1==pos2) return;
	for (size_t i=1;i<=line_count;i++)
		std::swap((*this)(pos1,i),(*this)(pos2,i));
}

template<class Type>
inline void Cmatrix<Type>::Lswap(const size_t pos1,const size_t pos2)
{
	if(pos1==0||pos2==0) throw bad_oper("Lswap出错:列标不能为0!");
	if(pos1>line_count||pos2>line_count) throw bad_oper("Lswap出错:列标越界!");
	if(pos1==pos2) return;
	matrix[pos1-1].swap(matrix[pos2-1]);
}

template<class Type>
inline void Cmatrix<Type>::Mscale(const Type factor,bool flag=true)
{
	if(flag)
	(*this)*=factor;
	else
	(*this)/=factor;
}

template<class Type>
inline void Cmatrix<Type>::Rscale(const size_t pos,const Type factor,bool flag=true)
{
	if(pos==0) throw bad_oper("Rscale出错:行标不能为0!");
	if(pos>row_count) throw bad_oper("Rscale出错:行标越界!");
	if(flag)
	{
		for (size_t i=1;i<=line_count;i++)
		(*this)(pos,i)*=factor;
	}
	else
	{
		for (size_t i=1;i<=line_count;i++)
		(*this)(pos,i)/=factor;
	}
}

template<class Type>
inline void Cmatrix<Type>::Lscale(const size_t pos,const Type factor,bool flag=true)
{
	if(pos==0) throw bad_oper("Lscale出错:列标不能为0!");
	if(pos>line_count) throw bad_oper("Lscale出错:列标越界!");
	if(flag)
	for (size_t i=1;i<=row_count;i++)
		(*this)(i,pos)*=factor;
	else
	for (size_t i=1;i<=row_count;i++)
		(*this)(i,pos)/=factor;
}

template<class Type>
inline void Cmatrix<Type>::Rsort(size_t beg,size_t end,Method flag=INC)
{
	if(end==0||beg==0) throw bad_oper("Rsort出错:行标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>row_count) throw bad_oper("Rsort出错:行标越界!");
	if(flag==INC)
	{
		for (size_t i=beg;i<=end;i++)
		{
			for (size_t j=line_count;j>=2;j--)
				for(size_t k=2;k<=j;k++)
					if((*this)(i,k-1)>(*this)(i,k)) std::swap((*this)(i,k-1),(*this)(i,k));
		}
	}
	else
	{
		for (size_t i=beg;i<=end;i++)
		{
			for (size_t j=line_count;j>=2;j--)
				for(size_t k=2;k<=j;k++)
					if((*this)(i,k-1)<(*this)(i,k)) std::swap((*this)(i,k-1),(*this)(i,k));
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::Lsort(size_t beg,size_t end,Method flag=INC)
{
	if(end==0||beg==0) throw bad_oper("Lsort出错:列标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>line_count) throw bad_oper("Lsort出错:列标越界!");
	if(flag==INC)
	{
		for (size_t i=beg;i<=end;i++)
		{
			for (size_t j=row_count;j>=2;j--)
				for(size_t k=2;k<=j;k++)
					if((*this)(k-1,i)>(*this)(k,i)) std::swap((*this)(k-1,i),(*this)(k,i));
		}
	}
	else
	{
		for (size_t i=beg;i<=end;i++)
			{
				for (size_t j=row_count;j>=2;j--)
					for(size_t k=2;k<=j;k++)
						if((*this)(k-1,i)<(*this)(k,i)) std::swap((*this)(k-1,i),(*this)(k,i));
			}
	}	
}

template<class Type>
inline void Cmatrix<Type>::Msort(Method flag1=LINE,Method flag2=INC)
{
	if(flag1==LINE)
	{
		if(flag2==INC)
		{
			for (size_t i=line_count*row_count;i>=2;i--)
				for(size_t j=2;j<=i;j++)
				{
					size_t line1=(j-1)/row_count+1,row1=(j-1)%row_count,line2=j/row_count+1,row2=j%row_count;
					if(row1==0)
					{
						row1=row_count;
						line1-=1;
					}
					if(row2==0)
					{
						row2=row_count;
						line2-=1;
					}
					if((*this)(row1,line1)>(*this)(row2,line2)) std::swap((*this)(row1,line1),(*this)(row2,line2));
				}
		}
		else
		{
			for (size_t i=line_count*row_count;i>=2;i--)
				for(size_t j=2;j<=i;j++)
				{
					size_t line1=(j-1)/row_count+1,row1=(j-1)%row_count,line2=j/row_count+1,row2=j%row_count;
					if(row1==0)
					{
						row1=row_count;
						line1-=1;
					}
					if(row2==0)
					{
						row2=row_count;
						line2-=1;
					}
					if((*this)(row1,line1)<(*this)(row2,line2)) std::swap((*this)(row1,line1),(*this)(row2,line2));
				}
		}
	}
	else
	{
		if(flag2==INC)
		{
			for (size_t i=line_count*row_count;i>=2;i--)
				for(size_t j=2;j<=i;j++)
				{
					size_t line1=(j-1)%line_count,row1=(j-1)/line_count+1,line2=j%line_count,row2=j/line_count+1;
					if(line1==0)
					{
						line1=line_count;
						row1-=1;
					}
					if(line2==0)
					{
						line2=line_count;
						row2-=1;
					}
					if((*this)(row1,line1)>(*this)(row2,line2)) std::swap((*this)(row1,line1),(*this)(row2,line2));
				}
		}
		else
		{
			for (size_t i=line_count*row_count;i>=2;i--)
				for(size_t j=2;j<=i;j++)
				{
					size_t line1=(j-1)%line_count,row1=(j-1)/line_count+1,line2=j%line_count,row2=j/line_count+1;
					if(line1==0)
					{
						line1=line_count;
						row1-=1;
					}
					if(line2==0)
					{
						line2=line_count;
						row2-=1;
					}
					if((*this)(row1,line1)<(*this)(row2,line2)) std::swap((*this)(row1,line1),(*this)(row2,line2));
				}
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::Rsorta(size_t beg,size_t end,Method flag=INC)
{
	if(end==0||beg==0) throw bad_oper("Rsorta出错:行标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>row_count) throw bad_oper("Rsorta出错:行标越界!");
	if(flag==INC)
	{
		for (size_t i=beg;i<=end;i++)
		{
			for (size_t j=line_count;j>=2;j--)
				for(size_t k=2;k<=j;k++)
					if(aufunc::abs((*this)(i,k-1))>aufunc::abs((*this)(i,k))) std::swap((*this)(i,k-1),(*this)(i,k));
		}
	}
	else
	{
		for (size_t i=beg;i<=end;i++)
		{
			for (size_t j=line_count;j>=2;j--)
				for(size_t k=2;k<=j;k++)
					if(aufunc::abs((*this)(i,k-1))<aufunc::abs((*this)(i,k))) std::swap((*this)(i,k-1),(*this)(i,k));
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::Lsorta(size_t beg,size_t end,Method flag=INC)
{
	if(end==0||beg==0) throw bad_oper("Lsorta出错:列标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(end>line_count) throw bad_oper("Lsorta出错:列标越界!");
	if(flag==INC)
	{
		for (size_t i=beg;i<=end;i++)
		{
			for (size_t j=row_count;j>=2;j--)
				for(size_t k=2;k<=j;k++)
					if(aufunc::abs((*this)(k-1,i))>aufunc::abs((*this)(k,i))) std::swap((*this)(k-1,i),(*this)(k,i));
		}
	}
	else
	{
		for (size_t i=beg;i<=end;i++)
			{
				for (size_t j=row_count;j>=2;j--)
					for(size_t k=2;k<=j;k++)
						if(aufunc::abs((*this)(k-1,i))<aufunc::abs((*this)(k,i))) std::swap((*this)(k-1,i),(*this)(k,i));
			}
	}	
}

template<class Type>
inline void Cmatrix<Type>::Msorta(Method flag1=LINE,Method flag2=INC)
{
	if(flag1==LINE)
	{
		if(flag2==INC)
		{
			for (size_t i=line_count*row_count;i>=2;i--)
				for(size_t j=2;j<=i;j++)
				{
					size_t line1=(j-1)/row_count+1,row1=(j-1)%row_count,line2=j/row_count+1,row2=j%row_count;
					if(row1==0)
					{
						row1=row_count;
						line1-=1;
					}
					if(row2==0)
					{
						row2=row_count;
						line2-=1;
					}
					if(aufunc::abs((*this)(row1,line1))>aufunc::abs((*this)(row2,line2))) std::swap((*this)(row1,line1),(*this)(row2,line2));
				}
		}
		else
		{
			for (size_t i=line_count*row_count;i>=2;i--)
				for(size_t j=2;j<=i;j++)
				{
					size_t line1=(j-1)/row_count+1,row1=(j-1)%row_count,line2=j/row_count+1,row2=j%row_count;
					if(row1==0)
					{
						row1=row_count;
						line1-=1;
					}
					if(row2==0)
					{
						row2=row_count;
						line2-=1;
					}
					if(aufunc::abs((*this)(row1,line1))<aufunc::abs((*this)(row2,line2))) std::swap((*this)(row1,line1),(*this)(row2,line2));
				}
		}
	}
	else
	{
		if(flag2==INC)
		{
			for (size_t i=line_count*row_count;i>=2;i--)
				for(size_t j=2;j<=i;j++)
				{
					size_t line1=(j-1)%line_count,row1=(j-1)/line_count+1,line2=j%line_count,row2=j/line_count+1;
					if(line1==0)
					{
						line1=line_count;
						row1-=1;
					}
					if(line2==0)
					{
						line2=line_count;
						row2-=1;
					}
					if(aufunc::abs((*this)(row1,line1))>aufunc::abs((*this)(row2,line2))) std::swap((*this)(row1,line1),(*this)(row2,line2));
				}
		}
		else
		{
			for (size_t i=line_count*row_count;i>=2;i--)
				for(size_t j=2;j<=i;j++)
				{
					size_t line1=(j-1)%line_count,row1=(j-1)/line_count+1,line2=j%line_count,row2=j/line_count+1;
					if(line1==0)
					{
						line1=line_count;
						row1-=1;
					}
					if(line2==0)
					{
						line2=line_count;
						row2-=1;
					}
					if(aufunc::abs((*this)(row1,line1))<aufunc::abs((*this)(row2,line2))) std::swap((*this)(row1,line1),(*this)(row2,line2));
				}
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::sort(const size_t *index,size_t num,Method flag1=INC,Method flag2=LINE)
{
	if(flag2==LINE)
	{
		if(flag1==INC)
		{
			for(size_t i=row_count;i>=2;--i)
			{
				for(size_t j=1;j<i;j++)
				{
					for(size_t k=0;k<num;k++)
					{
						if((*this)(j,index[k])>(*this)(j+1,index[k]))
						{
							Rswap(j,j+1);
							break;
						}
						else if((*this)(j,index[k])<(*this)(j+1,index[k]))
							break;
					}
				}
			}
		}
		else
		{
			for(size_t i=row_count;i>=2;--i)
			{
				for(size_t j=1;j<i;j++)
				{
					for(size_t k=0;k<num;k++)
					{
						if((*this)(j,index[k])<(*this)(j+1,index[k]))
						{
							Rswap(j,j+1);
							break;
						}
						else if((*this)(j,index[k])>(*this)(j+1,index[k]))
							break;
					}
				}
			}
		}
	}
	else
	{
		if(flag1==INC)
		{
			for(size_t i=line_count;i>=2;--i)
			{
				for(size_t j=1;j<i;j++)
				{
					for(size_t k=0;k<num;k++)
					{
						if((*this)(index[k],j)>(*this)(index[k],j+1))
						{
							Lswap(j,j+1);
							break;
						}
						else if((*this)(index[k],j)<(*this)(index[k],j+1))
							break;
					}
				}
			}
		}
		else
		{
			for(size_t i=line_count;i>=2;--i)
			{
				for(size_t j=1;j<i;j++)
				{
					for(size_t k=0;k<num;k++)
					{
						if((*this)(index[k],j)<(*this)(index[k],j+1))
						{
							Lswap(j,j+1);
							break;
						}
						else if((*this)(index[k],j)>(*this)(index[k],j+1))
							break;
					}
				}
			}
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::sorta(const size_t *index,size_t num,Method flag1=INC,Method flag2=LINE)
{
	if(flag2==LINE)
	{
		if(flag1==INC)
		{
			for(size_t i=row_count;i>=2;--i)
			{
				for(size_t j=1;j<i;j++)
				{
					for(size_t k=0;k<num;k++)
					{
						if(aufunc::abs((*this)(j,index[k]))>aufunc::abs((*this)(j+1,index[k])))
						{
							Rswap(j,j+1);
							break;
						}
						else if(aufunc::abs((*this)(j,index[k]))<aufunc::abs((*this)(j+1,index[k])))
							break;
					}
				}
			}
		}
		else
		{
			for(size_t i=row_count;i>=2;--i)
			{
				for(size_t j=1;j<i;j++)
				{
					for(size_t k=0;k<num;k++)
					{
						if(aufunc::abs((*this)(j,index[k]))<aufunc::abs((*this)(j+1,index[k])))
						{
							Rswap(j,j+1);
							break;
						}
						else if(aufunc::abs((*this)(j,index[k]))>aufunc::abs((*this)(j+1,index[k])))
							break;
					}
				}
			}
		}
	}
	else
	{
		if(flag1==INC)
		{
			for(size_t i=line_count;i>=2;--i)
			{
				for(size_t j=1;j<i;j++)
				{
					for(size_t k=0;k<num;k++)
					{
						if(aufunc::abs((*this)(index[k],j))>aufunc::abs((*this)(index[k],j+1)))
						{
							Lswap(j,j+1);
							break;
						}
						else if(aufunc::abs((*this)(index[k],j))<aufunc::abs((*this)(index[k],j+1)))
							break;
					}
				}
			}
		}
		else
		{
			for(size_t i=line_count;i>=2;--i)
			{
				for(size_t j=1;j<i;j++)
				{
					for(size_t k=0;k<num;k++)
					{
						if(aufunc::abs((*this)(index[k],j))<aufunc::abs((*this)(index[k],j+1)))
						{
							Lswap(j,j+1);
							break;
						}
						else if(aufunc::abs((*this)(index[k],j)>(*this)(index[k],j+1)))
							break;
					}
				}
			}
		}
	}
}


template<class Type>
inline void Cmatrix<Type>::_rand(Type v1,Type v2,Method flag) //用于非浮点数
{
	if(v1==v2) {*this=v1;return;}
	if(v1>v2) std::swap(v1,v2);
	std::srand(time(0));
	std::rand();
	if(flag==RANDSYM||flag==RANDHEMIT)
	{
		if(row_count!=line_count) throw bad_oper("rand出错:矩阵不为方阵!");
		for (size_t i=1;i<=row_count;i++)
		{
			for (size_t j=1;j<i;j++)
			{
				const Type value=std::rand();
				const Type ys=value%(v2-v1);
				const int cs=int((value-ys)/(v2-v1));
				if(ys!=0)
					(*this)(j,i)=(*this)(i,j)=ys+v1;
				else
				{
					if(cs%2==0) (*this)(j,i)=(*this)(i,j)=v1;
					else (*this)(j,i)=(*this)(i,j)=v2;
				}
			}
			const Type value=std::rand();
			const Type ys=value%(v2-v1);
			const Type cs=(value-ys)/(v2-v1);
			if(ys!=0)
				(*this)(i,i)=ys+v1;
			else
			{
				if(cs%2==0) (*this)(i,i)=v1;
				else (*this)(i,i)=v2;
			}
		}
	}
	else
	{
		for (size_t i=1;i<=row_count;++i)
			for (size_t j=1;j<=line_count;++j)
			{
				const Type value=std::rand();
				const Type ys=value%(v2-v1);
				const Type cs=(value-ys)/(v2-v1);
				if(ys!=0)
					(*this)(i,j)=ys+v1;
				else
				{
					if(cs%2==0) (*this)(i,j)=v1;
					else (*this)(i,j)=v2;
				}
			}
	}
}

template<class Type>
inline void Cmatrix<Type>::_randc(Type v1,Type v2,Method flag) //用于long,int,short的复数
{
	typedef Type::value_type type;
	std::srand(time(0));
	std::rand();
	type real1=v1.real(),real2=v2.real(),imag1=v1.imag(),imag2=v2.imag(),_real,_imag;
	if (real1>real2) std::swap(real1,real2);
	if (imag1>imag2) std::swap(imag1,imag2);
	if(flag==RANDSYM)
	{
		if(row_count!=line_count) throw bad_oper("rand出错:矩阵不为方阵!");
		for (size_t i=1;i<=row_count;i++)
		{
			for (size_t j=1;j<i;j++)
			{
				//求实部
				if(real1==real2) _real=real1;
				else
				{
					const type value=std::rand();
					const type ys=value%(real2-real1);
					const type cs=(value-ys)/(real2-real1);
					if(ys!=0)
						_real=ys+real1;
					else
					{
						if(cs%2==0) _real=real1;
						else _real=real2;
					}
				}
				//求虚部
				if(imag1==imag2) _imag=imag1;
				else
				{
					const type value=std::rand();
					const type ys=value%(imag2-imag1);
					const type cs=(value-ys)/(imag2-imag1);
					if(ys!=0)
						_imag=ys+imag1;
					else
					{
						if(cs%2==0) _imag=imag1;
						else _imag=imag2;
					}
				}
				(*this)(i,j)=(*this)(j,i)=Type(_real,_imag);
			}
			if(real1==real2) _real=real1;
			else
			{
				const type value=std::rand();
				const type ys=value%(real2-real1);
				const type cs=(value-ys)/(real2-real1);
				if(ys!=0)
					_real=ys+real1;
				else
				{
					if(cs%2==0) _real=real1;
					else _real=real2;
				}
			}
				//求虚部
			if(imag1==imag2) _imag=imag1;
			else
			{
				const type value=std::rand();
				const type ys=value%(imag2-imag1);
				const type cs=(value-ys)/(imag2-imag1);
				if(ys!=0)
					_imag=ys+imag1;
				else
				{
					if(cs%2==0) _imag=imag1;
					else _imag=imag2;
				}
			}
			(*this)(i,i)=Type(_real,_imag);
		}
	}
	else if(flag==RANDHEMIT)
	{
		if(row_count!=line_count) throw bad_oper("rand出错:矩阵不为方阵!");
		for (size_t i=1;i<=row_count;i++)
		{
			for (size_t j=1;j<i;j++)
			{
				//求实部
				if(real1==real2) _real=real1;
				else
				{
					const type value=std::rand();
					const type ys=value%(real2-real1);
					const type cs=(value-ys)/(real2-real1);
					if(ys!=0)
						_real=ys+real1;
					else
					{
						if(cs%2==0) _real=real1;
						else _real=real2;
					}
				}
				//求虚部
				if(imag1==imag2) _imag=imag1;
				else
				{
					const type value=std::rand();
					const type ys=value%(imag2-imag1);
					const type cs=(value-ys)/(imag2-imag1);
					if(ys!=0)
						_imag=ys+imag1;
					else
					{
						if(cs%2==0) _imag=imag1;
						else _imag=imag2;
					}
				}
				(*this)(i,j)=Type(_real,_imag);
				(*this)(j,i)=-(*this)(i,j);
			}
			if(real1==real2) _real=real1;
			else
			{
				const type value=std::rand();
				const type ys=value%(real2-real1);
				const type cs=(value-ys)/(real2-real1);
				if(ys!=0)
					_real=ys+real1;
				else
				{
					if(cs%2==0) _real=real1;
					else _real=real2;
				}
			}
			//求虚部
			if(imag1==imag2) _imag=imag1;
			else
			{
				const type value=std::rand();
				const type ys=value%(imag2-imag1);
				const type cs=(value-ys)/(imag2-imag1);
				if(ys!=0)
					_imag=ys+imag1;
				else
				{
					if(cs%2==0) _imag=imag1;
					else _imag=imag2;
				}
			}
			(*this)(i,i)=Type(_real,_imag);
		}
	}
	else
	{
		for (size_t i=1;i<=row_count;++i)
			for (size_t j=1;j<=line_count;++j)
			{
				//求实部
				if(real1==real2) _real=real1;
				else
				{
					const type value=std::rand();
					const type ys=value%(real2-real1);
					const type cs=(value-ys)/(real2-real1);
					if(ys!=0)
						_real=ys+real1;
					else
					{
						if(cs%2==0) _real=real1;
						else _real=real2;
					}
				}
				//求虚部
				if(imag1==imag2) _imag=imag1;
				else
				{
					const type value=std::rand();
					const type ys=value%(imag2-imag1);
					const type cs=(value-ys)/(imag2-imag1);
					if(ys!=0)
						_imag=ys+imag1;
					else
					{
						if(cs%2==0) _imag=imag1;
						else _imag=imag2;
					}
				}
				(*this)(i,j)=Type(_real,_imag);
			}
	}
}

template<class Type>
inline void Cmatrix<Type>::_randcu(Type v1,Type v2,Method flag) //用于unsigned的复数
{
	typedef Type::value_type type;
	std::srand(time(0));
	std::rand();
	type real1=v1.real(),real2=v2.real(),imag1=v1.imag(),imag2=v2.imag(),_real,_imag;
	if (real1>real2) std::swap(real1,real2);
	if (imag1>imag2) std::swap(imag1,imag2);
	if(flag==RANDSYM)
	{
		if(row_count!=line_count) throw bad_oper("rand出错:矩阵不为方阵!");
		for (size_t i=1;i<=row_count;i++)
		{
			for (size_t j=1;j<i;j++)
			{
				//求实部
				if(real1==real2) _real=real1;
				else
				{
					const type value=std::rand();
					const type ys=value%(real2-real1);
					const type cs=(value-ys)/(real2-real1);
					if(ys!=0)
						_real=ys+real1;
					else
					{
						if(cs%2==0) _real=real1;
						else _real=real2;
					}
				}
				//求虚部
				if(imag1==imag2) _imag=imag1;
				else
				{
					const type value=std::rand();
					const type ys=value%(imag2-imag1);
					const type cs=(value-ys)/(imag2-imag1);
					if(ys!=0)
						_imag=ys+imag1;
					else
					{
						if(cs%2==0) _imag=imag1;
						else _imag=imag2;
					}
				}
				(*this)(i,j)=(*this)(j,i)=Type(_real,_imag);
			}
			if(real1==real2) _real=real1;
			else
			{
				const type value=std::rand();
				const type ys=value%(real2-real1);
				const type cs=(value-ys)/(real2-real1);
				if(ys!=0)
					_real=ys+real1;
				else
				{
					if(cs%2==0) _real=real1;
					else _real=real2;
				}
			}
				//求虚部
			if(imag1==imag2) _imag=imag1;
			else
			{
				const type value=std::rand();
				const type ys=value%(imag2-imag1);
				const type cs=(value-ys)/(imag2-imag1);
				if(ys!=0)
					_imag=ys+imag1;
				else
				{
					if(cs%2==0) _imag=imag1;
					else _imag=imag2;
				}
			}
			(*this)(i,i)=Type(_real,_imag);
		}
	}
	else
	{
		for (size_t i=1;i<=row_count;++i)
			for (size_t j=1;j<=line_count;++j)
			{
				//求实部
				if(real1==real2) _real=real1;
				else
				{
					const type value=std::rand();
					const type ys=value%(real2-real1);
					const type cs=(value-ys)/(real2-real1);
					if(ys!=0)
						_real=ys+real1;
					else
					{
						if(cs%2==0) _real=real1;
						else _real=real2;
					}
				}
				//求虚部
				if(imag1==imag2) _imag=imag1;
				else
				{
					const type value=std::rand();
					const type ys=value%(imag2-imag1);
					const type cs=(value-ys)/(imag2-imag1);
					if(ys!=0)
						_imag=ys+imag1;
					else
					{
						if(cs%2==0) _imag=imag1;
						else _imag=imag2;
					}
				}
				(*this)(i,j)=Type(_real,_imag);
			}
	}
}

template<class Type>
inline void Cmatrix<Type>::rand(Type v1,Type v2,Method flag)
{
	srand(time(0));
	std::rand();
	if(flag==RANDSYM)
	{
		if(row_count!=line_count) throw bad_oper("rand出错:矩阵不为方阵!");
		for (size_t i=1;i<=row_count;i++)
		{
			for (size_t j=1;j<i;j++)
			{
				Type ratio=Type(std::rand())/Type(RAND_MAX);
				(*this)(i,j)=(v2-v1);
				(*this)(i,j)*=ratio;
				(*this)(i,j)+=v1;
				(*this)(j,i)=(*this)(i,j);
			}
			Type ratio=Type(std::rand())/Type(RAND_MAX);
			(*this)(i,i)=(v2-v1);
			(*this)(i,i)*=ratio;
			(*this)(i,i)+=v1;
		}
	}
	else if(flag==RANDHEMIT)
	{
		if(row_count!=line_count) throw bad_oper("rand出错:矩阵不为方阵!");
		for (size_t i=1;i<=row_count;i++)
		{
			for (size_t j=1;j<i;j++)
			{
				Type ratio=Type(std::rand())/Type(RAND_MAX);
				(*this)(i,j)=(v2-v1);
				(*this)(i,j)*=ratio;
				(*this)(i,j)+=v1;
				(*this)(j,i)=aufunc::conj((*this)(i,j));
			}
			Type ratio=Type(std::rand())/Type(RAND_MAX);
			(*this)(i,i)=(v2-v1);
			(*this)(i,i)*=ratio;
			(*this)(i,i)+=v1;
		}
	}
	else
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
		{
			Type ratio=Type(std::rand())/Type(RAND_MAX);
			(*this)(i,j)=(v2-v1);
			(*this)(i,j)*=ratio;
			(*this)(i,j)+=v1;
		}
}

template<>
inline void Cmatrix<int>::rand(int v1,int v2,Method flag)
{_rand(v1,v2,flag);}

template<>
inline void Cmatrix<short>::rand(short v1,short v2,Method flag)
{_rand(v1,v2,flag);}

template<>
inline void Cmatrix<long>::rand(long v1,long v2,Method flag)
{_rand(v1,v2,flag);}

template<>
inline void Cmatrix<unsigned int>::rand(unsigned int v1,unsigned int v2,Method flag)
{_rand(v1,v2,flag);}

template<>
inline void Cmatrix<unsigned short>::rand(unsigned short v1,unsigned short v2,Method flag)
{_rand(v1,v2,flag);}

template<>
inline void Cmatrix<unsigned long>::rand(unsigned long v1,unsigned long v2,Method flag)
{_rand(v1,v2,flag);}

template<>
inline void Cmatrix<std::complex<int> >::rand(std::complex<int> v1,std::complex<int> v2,Method flag)
{
	_randc(v1,v2,flag);
}

template<>
inline void Cmatrix<std::complex<short> >::rand(std::complex<short> v1,std::complex<short> v2,Method flag)
{
	_randc(v1,v2,flag);
}

template<>
inline void Cmatrix<std::complex<long> >::rand(std::complex<long> v1,std::complex<long> v2,Method flag)
{
	_randc(v1,v2,flag);
}

template<>
inline void Cmatrix<std::complex<unsigned long> >::rand(std::complex<unsigned long> v1,std::complex<unsigned long> v2,Method flag)
{
	_randcu(v1,v2,flag);
}

template<>
inline void Cmatrix<std::complex<unsigned short> >::rand(std::complex<unsigned short> v1,std::complex<unsigned short> v2,Method flag)
{
	_randcu(v1,v2,flag);
}

template<>
inline void Cmatrix<std::complex<unsigned int> >::rand(std::complex<unsigned int> v1,std::complex<unsigned int> v2,Method flag)
{
	_randcu(v1,v2,flag);
}


template<class Type>
inline void Cmatrix<Type>::clear()
{
	matrix.clear();
	line_count=row_count=0;
}

template<class Type>
inline void Cmatrix<Type>::Lclear()
{
	matrix.clear();
	line_count=0;
}

template<class Type>
inline void Cmatrix<Type>::Rclear()
{
	for(size_t i=0;i<line_count;++i)
		matrix[i].clear();
	row_count=0;
}

//赋值操作************************************************
template<class Type>
template<class U>
inline void Cmatrix<Type>::assign(const U *ptr,size_type num,Method flag=LINE)
{
	if (flag==LINE)
	{
		for (size_type i=0;i<line_count;i++)
			for (size_type j=0;j<row_count;j++)
				if(num-->0)
				matrix[i][j]=*ptr++;
				else
				return;
	}
	else
	{
		for (size_type i=0;i<row_count;i++)
			for (size_type j=0;j<line_count;j++)
				if(num--)
				matrix[j][i]=*ptr++;
				else
				return;
	}
};

template<class Type>
template <class InputIterator>
inline void Cmatrix<Type>::assign(InputIterator beg,InputIterator end,Method flag=LINE)
{
	if (flag==LINE)
	{
		for (size_type i=0;i<line_count;i++)
			for (size_type j=0;j<row_count;j++)
				if(beg!=end)
				matrix[i][j]=*beg++;
				else return;
	}
	else
	{
		for (size_type i=0;i<row_count;i++)
			for (size_type j=0;j<line_count;j++)
				if(beg!=end)
				matrix[j][i]=*beg++;
				else return;
	}
};

template<class Type>
template <class InputIterator>
inline void Cmatrix<Type>::assign(InputIterator beg,Method flag=LINE)
{
	if (flag==LINE)
	{
		for (size_type i=0;i<line_count;i++)
			for (size_type j=0;j<row_count;j++)
				matrix[i][j]=*beg++;
	}
	else
	{
		for (size_type i=0;i<row_count;i++)
			for (size_type j=0;j<line_count;j++)
				matrix[j][i]=*beg++;
	}
};

template<class Type>
template<class U>
inline void Cmatrix<Type>::assign(size_t row1,size_t line1,size_t row2,size_t line2, const Cmatrix<U> &orig,const size_t row3,const size_t line3)
{
	if(row1==0||row2==0||row3==0)	throw bad_oper("assign出错:行标不能为0!");
	if(line1==0||line2==0||line3==0)	throw bad_oper("assign出错:列标不能为0!");
	if(row1>row2) std::swap(row1,row2);
	if(line1>line2) std::swap(line1,line2);
	if(row2>row_count) throw bad_oper("assign出错:行标越界!");
	if(line2>line_count) throw bad_oper("assign出错:列标越界!");
	const size_t row4=row3+row2-row1,line4=line3+line2-line1;
	if(row4>orig.row()) throw bad_oper("assign出错:行标越界!");
	if(line4>orig.line()) throw bad_oper("assign出错:列标越界!");
	for(size_t i=row1;i<=row2;i++)
		for(size_t j=line1;j<=line2;j++)
			(*this)(i,j)=orig(i-row1+row3,j-line1+line3);
}

template <class Type> template <class U>
inline Cmatrix<Type> & Cmatrix<Type>::operator =(const Cmatrix<U> &origin)
{
	if(this==&origin) return *this;
	resize(origin.row(),origin.line());
	for (size_type row=1;row<=row_count;row++)
		for (size_type line=1;line<=line_count;line++)
			matrix[line-1][row-1]=origin(row,line);
	return *this;
};

template <class Type> 
inline Cmatrix<Type> & Cmatrix<Type>::operator =(const Type a)
{
	if(IsEmpty()) throw bad_oper("operator+=出错:矩阵为空!");
	for (size_type line=0;line<line_count;++line)
		for (size_type row=0;row<row_count;++row)
			matrix[line][row]=a;
	return *this;
};

//数学运算符的重载1************************************************
template <class	Type> template <class U>
inline Cmatrix<Type>& Cmatrix<Type>::operator +=(const Cmatrix<U> &a)
{
	if(IsEmpty()) throw bad_oper("operator+=出错:矩阵为空!");
	if (line_count!=a.line()||row_count!=a.row())
		throw bad_oper("operator+=出错:两矩阵行列数不等!");
	for (size_type line=0;line<line_count;++line)
		for (size_type row=0;row<row_count;++row)
			matrix[line][row]+=a(row+1,line+1);
	return *this;
};

template <class	Type> template <class U>
inline Cmatrix<Type>& Cmatrix<Type>::operator -=(const Cmatrix<U> &a)
{
	if(IsEmpty()) throw bad_oper("operator-=出错:矩阵为空!");
	if (line_count!=a.line()||row_count!=a.row())
		throw bad_oper("operator-=出错:两矩阵行列数不等!");
	for (size_type line=0;line<line_count;++line)
		for (size_type row=0;row<row_count;++row)
			matrix[line][row]-=a(row+1,line+1);
	return *this;
};

template <class	Type> template <class U>
inline Cmatrix<Type>& Cmatrix<Type>::operator *=(const Cmatrix<U> &a)
{
	if(IsEmpty()||a.IsEmpty()) throw bad_oper("operator*=出错:矩阵为空!");
	if (line_count!=a.row())
		throw bad_oper("operator*=出错:两矩阵无法相乘!");
	std::vector<Type> tmp(a.line());
	if (a.line()>line_count)
	{
		matrix.resize(a.line());
		for (size_type line=line_count;line<a.line();matrix[line++].resize(row_count));
	}
	for (size_type row=1;row<=row_count;row++)
	{
		tmp.assign(a.line(),Type());
		for (size_type line=1;line<=a.line();line++)
		{
			for (size_type i=1;i<=line_count;i++)
				tmp[line-1]+=(*this)(row,i)*a(i,line);
		}
		for (size_type line=0;line<a.line();line++)
			matrix[line][row-1]=tmp[line];
	}
	Lresize(a.line());
	return *this;
};

template <class	Type> template <class U>
inline Cmatrix<Type>& Cmatrix<Type>::operator /=(const Cmatrix<U> &a)
{
	return (*this)*=a.inver();
};

template <class	Type> template <class U>
inline Cmatrix<Type> Cmatrix<Type>::operator +(const Cmatrix<U> &a) const
{
	if(IsEmpty()) throw bad_oper("operator+出错:矩阵为空!");
	if (line_count!=a.line()||row_count!=a.row())
		throw bad_oper("operator+出错:两矩阵行列数不等!");
	Cmatrix<Type> temp(*this);
	for (size_type line=0;line<line_count;++line)
		for (size_type row=0;row<row_count;++row)
			temp.matrix[line][row]+=a(row+1,line+1);
	return temp;
};

template <class	Type> template <class U>
inline Cmatrix<Type> Cmatrix<Type>::operator -(const Cmatrix<U> &a) const
{
	if(IsEmpty()) throw bad_oper("operator-出错:矩阵为空!");
	if (line_count!=a.line()||row_count!=a.row())
		throw bad_oper("operator-出错:两矩阵行列数不等!");
	Cmatrix<Type> temp(*this);
	for (size_type line=0;line<line_count;++line)
		for (size_type row=0;row<row_count;++row)
			temp.matrix[line][row]-=a(row+1,line+1);
	return temp;
};

template <class	Type> template <class U>
inline Cmatrix<Type> Cmatrix<Type>::operator *(const Cmatrix<U> &a) const
{
	if(IsEmpty()||a.IsEmpty()) throw bad_oper("operator*出错:矩阵为空!");
	if (line_count!=a.row())
		throw bad_oper("operator*出错:两矩阵无法相乘!");
	Cmatrix<Type> temp(row_count,a.line());
	for (size_type row=1;row<=temp.row_count;row++)
		for (size_type line=1;line<=temp.line_count;line++)
		{
			for (size_type i=1;i<=line_count;i++)
				temp(row,line)+=matrix[i-1][row-1]*a(i,line);
		}
	return temp;
};

template <class	Type> template <class U>
inline Cmatrix<Type> Cmatrix<Type>::operator /(const Cmatrix<U> &a) const
{
		return a.inver()*=(*this);
};

//数学运算符的重载2************************************************
template <class	Type> 
inline Cmatrix<Type>& Cmatrix<Type>::operator +=(const Type a)
{
	if(IsEmpty()) throw bad_oper("operator+=出错:矩阵为空!");
	for (size_type line=0;line<line_count;++line)
		for (size_type row=0;row<row_count;++row)
			matrix[line][row]+=a;
	return *this;
};

template <class	Type> 
inline Cmatrix<Type>& Cmatrix<Type>::operator -=(const Type a)
{
	if(IsEmpty()) throw bad_oper("operator-=出错:矩阵为空!");
	for (size_type line=0;line<line_count;++line)
		for (size_type row=0;row<row_count;++row)
			matrix[line][row]-=a;
	return *this;
};

template <class	Type> 
inline Cmatrix<Type>& Cmatrix<Type>::operator *=(const Type a)
{
	if(IsEmpty()) throw bad_oper("operator*=出错:矩阵为空!");
	for (size_type line=0;line<line_count;++line)
		for (size_type row=0;row<row_count;++row)
			matrix[line][row]*=a;
	return *this;
};

template <class	Type> 
inline Cmatrix<Type>& Cmatrix<Type>::operator /=(const Type a)
{
	if(IsEmpty()) throw bad_oper("operator/=出错:矩阵为空!");
	for (size_type line=0;line<line_count;++line)
		for (size_type row=0;row<row_count;++row)
			matrix[line][row]/=a;
	return *this;
};

template <class	Type> 
inline Cmatrix<Type> Cmatrix<Type>::operator +(const Type a) const
{
	if(IsEmpty()) throw bad_oper("operator+出错:矩阵为空!");
	Cmatrix<Type> temp(*this);
	for (size_type line=0;line<line_count;++line)
		for (size_type row=0;row<row_count;++row)
			temp.matrix[line][row]+=a;
	return temp;
};

template <class	Type> 
inline Cmatrix<Type> Cmatrix<Type>::operator -(const Type a) const
{
	if(IsEmpty()) throw bad_oper("operator-出错:矩阵为空!");
	Cmatrix<Type> temp(*this);
	for (size_type line=0;line<line_count;++line)
		for (size_type row=0;row<row_count;++row)
			temp.matrix[line][row]-=a;
	return temp;
};

template <class	Type> 
inline Cmatrix<Type> Cmatrix<Type>::operator *(const Type a) const
{
	if(IsEmpty()) throw bad_oper("operator*出错:矩阵为空!");
	Cmatrix<Type> temp(*this);
	for (size_type row=1;row<=temp.row_count;row++)
		for (size_type line=1;line<=temp.line_count;line++)
				temp(row,line)*=a;
	return temp;
};

template <class	Type> 
inline Cmatrix<Type> Cmatrix<Type>::operator /(const Type a) const
{
	if(IsEmpty()) throw bad_oper("operator/出错:矩阵为空!");
	Cmatrix<Type> temp(*this);
	for (size_type row=1;row<=temp.row_count;row++)
		for (size_type line=1;line<=temp.line_count;line++)
				temp(row,line)/=a;
	return temp;
};

//输入输出
template<class Type, class elem, class traits>
inline std::basic_istream<elem, traits>&operator>>(std::basic_istream<elem, traits>& is,  Cmatrix<Type>&origin)
{
	if(origin.IsEmpty()) return is;
	for(size_t row=1;row<=origin.row();row++)
		for (size_t line=1;line<=origin.line();line++)
			is>>origin(row,line);
	return is;
};

template<class Type, class elem, class traits>
inline std::basic_ostream<elem, traits>&operator<< (std::basic_ostream<elem, traits>& os, const Cmatrix<Type>&origin)
{
	if(origin.IsEmpty()) return os;
	os<<origin(1,1);
	for (size_t line=2;line<=origin.line();line++)
		os<<' '<<origin(1,line);
	for(size_t row=2;row<=origin.row();row++)
	{
		os<<std::endl;
		os<<origin(row,1);
		for (size_t line=2;line<=origin.line();line++)
			os<<' '<<origin(row,line);
	}
	return os;
};


//变动形矩阵运算************************************************

//加
template<class Type>
template<class U>
inline void Cmatrix<Type>::add(size_t row1,size_t line1,size_t row2,size_t line2,const Cmatrix<U>&orig,const size_t row3,const size_t line3)
{
	if(row1==0||row2==0||row3==0)	throw bad_oper("add出错:行标不能为0!");
	if(line1==0||line2==0||line3==0)	throw bad_oper("add出错:列标不能为0!");
	if(row1>row2) std::swap(row1,row2);
	if(line1>line2) std::swap(line1,line2);
	if(row2>row_count) throw bad_oper("add出错:行标越界!");
	if(line2>line_count) throw bad_oper("add出错:列标越界!");
	const size_t row4=row3+row2-row1,line4=line3+line2-line1;
	if(row4>orig.row()) throw bad_oper("add出错:行标越界!");
	if(line4>orig.line()) throw bad_oper("add出错:列标越界!");
	for(size_t i=row1;i<=row2;i++)
		for(size_t j=line1;j<=line2;j++)
			(*this)(i,j)+=orig(i-row1+row3,j-line1+line3);
}

//减
template<class Type>
template<class U>
inline void Cmatrix<Type>::sub(size_t row1,size_t line1,size_t row2,size_t line2,const Cmatrix<U>&orig,const size_t row3=1,const size_t line3=1)
{
	if(row1==0||row2==0||row3==0)	throw bad_oper("sub出错:行标不能为0!");
	if(line1==0||line2==0||line3==0)	throw bad_oper("sub出错:列标不能为0!");
	if(row1>row2) std::swap(row1,row2);
	if(line1>line2) std::swap(line1,line2);
	if(row2>row_count) throw bad_oper("sub出错:行标越界!");
	if(line2>line_count) throw bad_oper("sub出错:列标越界!");
	const size_t row4=row3+row2-row1,line4=line3+line2-line1;
	if(row4>orig.row()) throw bad_oper("sub出错:行标越界!");
	if(line4>orig.line()) throw bad_oper("sub出错:列标越界!");
	for(size_t i=row1;i<=row2;i++)
		for(size_t j=line1;j<=line2;j++)
			(*this)(i,j)-=orig(i-row1+row3,j-line1+line3);
}

//点乘
template<class Type>
template<class U>
inline bool Cmatrix<Type>::Pmutiply(const Cmatrix<U>&R)
{
	if(IsEmpty()||R.IsEmpty()) return false;
	if (row_count!=R.row()||line_count!=R.line())
		return false;
		for (size_t i=1;i<=row_count;++i)
			for (size_t j=1;j<=line_count;++j)
				(*this)(i,j)*=R(i,j);
	return true;
}

//左叉乘
template<class Type>
template<class U>
inline void Cmatrix<Type>::LCmutiply(const Cmatrix<U>&R)
{
	if(IsEmpty()||R.IsEmpty()) throw bad_oper("LCmutiply出错:矩阵为空!");
	const size_t row=row_count,line=line_count;
	resize(row_count*R.row(),line_count*R.line());
	for (size_t i=row;i>=1;--i)
		for (size_t j=line;j>=1;--j)
		{
			for (size_t ii=i*R.row();ii>=(i-1)*R.row()+1;--ii)
				for (size_t jj=j*R.line();jj>=(j-1)*R.line()+1;--jj)
                (*this)(ii,jj)=(*this)(i,j)*R((ii-(i-1)*R.row()),(jj-(j-1)*R.line()));
		};
}

//右叉乘
template<class Type>
template<class U>
inline void Cmatrix<Type>::RCmutiply(const Cmatrix<U>&L)
{
	if(IsEmpty()||L.IsEmpty()) throw bad_oper("RCmutiply出错:矩阵为空!");
	const size_t row=row_count,line=line_count;
	resize(row_count*L.row(),line_count*L.line());
	for (size_t i=L.row();i>=1;--i)
		for (size_t j=L.line();j>=1;--j)
		{
			for (size_t ii=i*row;ii>=(i-1)*row+1;--ii)
				for (size_t jj=j*line;jj>=(j-1)*line+1;--jj)
                (*this)(ii,jj)=L(i,j)*(*this)((ii-(i-1)*row),(jj-(j-1)*line));
		}
}

//左乘
template<typename Type>
template<class U>
inline bool Cmatrix<Type>::Lmutiply(const Cmatrix<U> & orig)
{
	*this*=orig;
	return true;
}

//右乘
template<typename Type>
template<class U>
inline bool Cmatrix<Type>::Rmutiply(const Cmatrix<U>& orig)
{
	if(IsEmpty()||orig.IsEmpty()) return false;
	if (orig.line()!=row_count)	return false;
	std::vector<Type> tmp(orig.row());
	if (orig.row()>row_count) 
	{
		Rresize(orig.row());
	}
	for (size_type line=1;line<=line_count;line++)
	{
		tmp.assign(orig.row(),Type());
		for (size_type row=1;row<=orig.row();row++)
		{
			for (size_type i=1;i<=orig.line();i++)
				tmp[row-1]+=orig(row,i)*(*this)(i,line);
		}
		for (size_type row=0;row<orig.row();row++)
			matrix[line-1][row]=tmp[row];
	}
	Rresize(orig.row());
	return true;
}

//乘
template<typename Type>
template<class U,class V>
inline bool Cmatrix<Type>::mutiply(const Cmatrix<U>& L,const Cmatrix<V>& R)
{
	if(L.IsEmpty()||R.IsEmpty()) return false;
	if (L.line()!=R.row())	return false;
	resize(L.row(),R.line());
	for(size_t i=1;i<=row_count;i++)
		for(size_t j=1;j<=line_count;j++)
		{
			(*this)(i,j)=0;
			for(size_t k=1;k<=L.line();k++)
				(*this)(i,j)+=L(i,k)*R(k,j);
		}
	return true;
}

//单位化
template<class Type>
inline void Cmatrix<Type>::united(size_type num)
{
	resize(num,num);
	for (size_type row=1;row<=num;row++) 
		for (size_type line=1;line<=num;line++)
		{
			if(row==line)(*this)(row,line)=1;
			else (*this)(row,line)=0;
		}
}

//转置
template<class Type>
inline void Cmatrix<Type>::transed()
{
	if(IsEmpty()) throw bad_oper("transed出错:矩阵为空!");
	size_type tmp;
	tmp=(row_count<line_count?row_count:line_count);
	for(size_type row=1;row<=tmp;row++)
		for (size_type line=row+1;line<=tmp;line++)
			std::swap((*this)(row,line),(*this)(line,row));
	const int dist=line_count-row_count;
	if (dist>0)
	{
		for (size_type line=0;line<tmp;line++)
			for(size_type i=1;i<=dist;i++)
				matrix[line].push_back((*this)(line+1,tmp+i));
	}
	if (dist<0)
	{
		matrix.resize(row_count);
		for (size_type line=tmp;line<row_count;line++)
			for(size_type row=1;row<=tmp;row++)
				matrix[line].push_back((*this)(line+1,row));
	}
	resize(line_count,row_count);
};

//共轭转置
template<class Type>
inline void Cmatrix<Type>::hermited()
{
	if(IsEmpty()) throw bad_oper("hermited出错:矩阵为空!");;
	size_type tmp;
	tmp=(row_count<line_count?row_count:line_count);
	for(size_type row=1;row<=tmp;row++)
		for (size_type line=row+1;line<=tmp;line++)
			std::swap((*this)(row,line),(*this)(line,row));
	const int dist=line_count-row_count;
	if (dist>0)
	{
		for (size_type line=0;line<tmp;line++)
			for(size_type i=1;i<=dist;i++)
				matrix[line].push_back((*this)(line+1,tmp+i));
	}
	if (dist<0)
	{
		matrix.resize(row_count);
		for (size_type line=tmp;line<row_count;line++)
			for(size_type row=1;row<=tmp;row++)
				matrix[line].push_back((*this)(line+1,row));
	}
	resize(line_count,row_count);

	for (size_t i=0;i<line_count;++i)
		for (size_t j=0;j<row_count;++j)
			matrix[i][j]=aufunc::conj(matrix[i][j]);
};


//逆矩阵化
template<class Type>
inline bool Cmatrix<Type>::invered()
{
	if(IsEmpty()) return false;
	if(line_count!=row_count) return false;
	Cmatrix<Type> tmp(row_count,line_count*2);
	tmp=*this;
	tmp.combine(Cmatrix<Type>(UNIT,line_count));
	std::pair<Type,size_t> vmax_pos;
	size_t pos=1; //记录搜寻主元的起始行
	for (size_t i=1;i<=line_count;i++)
	{
		vmax_pos=tmp.Lmaxa(i,pos,tmp.row());
		if (aufunc::abs(vmax_pos.first)<precision) return false;
		tmp.Rswap(pos,vmax_pos.second);
		tmp.Rscale(pos,tmp(pos,i),false);
		for(size_t j=1;j<=tmp.row();j++)
		{
			if(j==pos) continue;
			tmp.Radd(pos,j,j,-tmp(j,i));
		}
		++pos;
	}
	assign(tmp.begin()+row_count*row_count,tmp.end());
	return true;
};

//翻转
template<class Type>
inline void Cmatrix<Type>::revered(Method flag=LINE)
{
	if(IsEmpty()) throw bad_oper("revered出错:矩阵为空!");;
	if(flag==LINE)
	{
		for (size_t i=1;i<=line_count;i++)
			for (size_t j=1;j<=row_count/2;j++)
				std::swap((*this)(j,i),(*this)(row_count-j+1,i));
	}
	else
	{
		for (size_t i=1;i<=line_count/2;i++)
		std::swap(matrix[i-1],matrix[line_count-i]);
	}
};

template<class Type>
inline void Cmatrix<Type>::stded()
{
	if(IsEmpty()) throw bad_oper("stded出错:矩阵为空!");;
	std::pair<Type,size_t> vmax_pos;
	size_t pos=1; //记录搜寻主元的其实行
	for (size_t i=1;i<=line_count;i++)
	{
		if(pos>row_count) return ;
		vmax_pos=Lmaxa(i,pos,row_count);
		if (aufunc::abs(vmax_pos.first)<precision) continue;
		Rswap(pos,vmax_pos.second);
		Rscale(pos,(*this)(pos,i),false);
		for(size_t j=pos+1;j<=row();j++)
			Radd(pos,j,j,-(*this)(j,i));
		for(size_t j=1;j<pos;j++)
			Radd(pos,j,j,-(*this)(j,i));
		++pos;
	}
}

template<class Type>
inline void Cmatrix<Type>::gaussed() 
{
	if(IsEmpty()) throw bad_oper("gaussed出错:矩阵为空!");;
	std::pair<Type,size_t> vmax_pos;
	size_t pos=1; //记录搜寻主元的其实行
	for (size_t i=1;i<=line_count;i++)
	{
		if(pos>row_count) return ;
		vmax_pos=Lmaxa(i,pos,row_count);
		if (aufunc::abs(vmax_pos.first)<precision) continue;
		Rswap(pos,vmax_pos.second);
		for(size_t j=pos+1;j<=row();j++)
			Radd(pos,j,j,-(*this)(j,i)/(*this)(pos,i));
		++pos;
	}
}

template<class Type>
inline void Cmatrix<Type>::zeroed()
{
	for(size_t i=1;i<=row_count;i++)
		for(size_t j=1;j<=line_count;j++)
			if(aufunc::abs((*this)(i,j))<precision) (*this)(i,j)=0;
}

template<class Type>
inline bool Cmatrix<Type>::orthorized(Method flag=LINE)
{
	if(IsEmpty()) return false;
	long double value=0;
	if(flag==LINE)
	{
		for(size_t i=1;i<=line_count;i++)
		{
			value=0;
			Cmatrix<Type> tmp(*this,i,i);
			for(size_t j=1;j<i;j++)
			{
				Type value1=0,value2=0;
				for(size_t k=1;k<=row_count;k++)
				{
					value1+=(*this)(k,j)*aufunc::conj(tmp(k,1));
					value2+=(*this)(k,j)*aufunc::conj((*this)(k,j));
				}
				for(size_t k=1;k<=row_count;k++)
					(*this)(k,i)-=value1/value2*(*this)(k,j);
			}
			for(size_t j=1;j<=row_count;j++)
				value+=aufunc::abs((*this)(j,i)*(*this)(j,i));
			if(value<precision) return false;
		}
		vec_united();
		return true;
	}
	else
	{
		for(size_t i=1;i<=row_count;i++)
		{
			value=0;
			Cmatrix<Type> tmp(*this,i,i,ROW);
			for(size_t j=1;j<i;j++)
			{
				Type value1=0,value2=0;
				for(size_t k=1;k<=line_count;k++)
				{
					value1+=(*this)(j,k)*aufunc::conj(tmp(1,k));
					value2+=(*this)(j,k)*aufunc::conj((*this)(j,k));
				}
				for(size_t k=1;k<=line_count;k++)
					(*this)(i,k)-=value1/value2*(*this)(j,k);
			}
			for(size_t j=1;j<=line_count;j++)
				value+=aufunc::abs((*this)(i,j)*(*this)(i,j));
			if(value<precision) return false;
		}
		vec_united(ROW);
		return true;
	}	
}

template<class Type>
inline void Cmatrix<Type>::vec_united(Method flag= LINE)
{
	if(IsEmpty()) throw bad_oper("vec_united出错:矩阵为空!");;
	if(flag==LINE)
	{
		for(size_t i=1;i<=line_count;i++)
		{
			long double value=0;
			for(size_t j=1;j<=row_count;j++)
				value+=aufunc::abs((*this)(j,i)*(*this)(j,i));
			for(size_t j=1;j<=row_count;j++)
				(*this)(j,i)/=std::pow(value,long double(0.5));
		}
	}
	else
	{
		for(size_t i=1;i<=row_count;i++)
		{
			long double value=0;
			for(size_t j=1;j<=line_count;j++)
				value+=aufunc::abs((*this)(i,j)*(*this)(i,j));
			for(size_t j=1;j<=line_count;j++)
				(*this)(i,j)/=std::pow(value,long double(0.5));
		}
	}
}

template<class Type>
inline void Cmatrix<Type>::vec_normalized(Method flag= LINE)
{
	if(IsEmpty()) throw bad_oper("vec_normalized出错:矩阵为空!");
	if(flag==LINE)
	{
		for(size_t i=1;i<=line_count;i++)
		{
			Type value=Lmaxa(i,1,row_count).first;
			for(size_t j=1;j<=row_count;j++)
				(*this)(j,i)/=value;
		}
	}
	else
	{
		for(size_t i=1;i<=row_count;i++)
		{
			Type value=Rmaxa(i,1,line_count).first;
			for(size_t j=1;j<=line_count;j++)
				(*this)(i,j)/=value;
		}
	}
}


template<class Type>
inline bool Cmatrix<Type>::eig_jacob (Cmatrix<Type>& eig,Cmatrix<Type>&vec,Method flag=INC)
{
	if(!issym()) return false;
	eig.resize(row_count,1);
	vec.united(row_count);
	long double value=0;
	const long double n=(row_count*line_count-row_count)/2.0;
	value=0;
	for(size_t i=1;i<=row_count;i++)
		for(size_t j=1;j<i;j++)
			value+=aufunc::abs((*this)(i,j)*(*this)(i,j));
	do
	{
		//计算非对角元素的平方和
		value/=n;
		size_t lpos=0,rpos=0;
		for(rpos=1;rpos<=row_count;rpos++)
		{
			for(lpos=rpos+1;lpos<=line_count;lpos++)
			{
				if(aufunc::abs((*this)(rpos,lpos))>value)
				{
					long double y=aufunc::abs((*this)(rpos,rpos)-(*this)(lpos,lpos));
					long double x=2.0*(*this)(rpos,lpos);
					if((*this)(rpos,rpos)<(*this)(lpos,lpos)) x*=(-1.0);
					long double cos2a=y/std::sqrt(x*x+y*y);
					long double cosa=std::sqrt((1+cos2a)/2);
					long double sin2a=x/std::sqrt(x*x+y*y);
					long double sina=sin2a/(2*cosa);
					const long double pq=(*this)(rpos,lpos),pp=(*this)(rpos,rpos),qq=(*this)(lpos,lpos);
					for(size_t j=1;j<=line_count;j++)
					{
						long double qj=(*this)(lpos,j),pj=(*this)(rpos,j);
						(*this)(j,rpos)=(*this)(rpos,j)=pj*cosa+qj*sina;
						(*this)(lpos,j)=(*this)(j,lpos)=(0-pj)*sina+qj*cosa;
					}
					(*this)(rpos,rpos)=pp*cosa*cosa+2*pq*sina*cosa+qq*sina*sina;
					(*this)(lpos,lpos)=pp*sina*sina-2*pq*sina*cosa+qq*cosa*cosa;
					(*this)(rpos,lpos)=(*this)(lpos,rpos)=0.5*(qq-pp)*sin2a+pq*cos2a;
					Cmatrix<Type> R(UNIT,row_count);
					R(rpos,rpos)=R(lpos,lpos)=cosa;
					R(rpos,lpos)=-sina;R(lpos,rpos)=sina;
					vec*=R;
				}
			}
		}
	}while(value>precision);
	for(size_t i=1;i<=row_count;i++)
		eig(i,1)=(*this)(i,i);
	if(flag==INC)
	{
		for(size_t i=1;i<=row_count;i++)
		{
			std::pair<Type,size_t> val_pos;
			val_pos=eig.Lmin(1,i,eig.row());
			eig.Rswap(val_pos.second,i);
			vec.Lswap(val_pos.second,i);
		}
	}
	else if(flag==DEC)
	{
		for(size_t i=1;i<=row_count;i++)
		{
			std::pair<Type,size_t> val_pos;
			val_pos=eig.Lmax(1,i,eig.row());
			eig.Rswap(val_pos.second,i);
			vec.Lswap(val_pos.second,i);
		}
	}
	return true;
}

template<class Type>
inline bool Cmatrix<Type>::eig_subspace(Cmatrix<Type>& eig,Cmatrix<Type>& vec,const size_t num=1)
{
	if(!issym()) return false;
	if(num>line_count) return false;
	eig.resize(num,1);
	vec.resize(row_count,num);
	Cmatrix<Type> C(row_count,num),B(num,num);
	for(size_t i=1;i<=num;i++)
		vec(i,i)=1;
	vec.orthorized();
	Cmatrix<Type> _eig(num);
	do
	{
		eig=_eig;
		C=vec;
		C.Rmutiply(*this);
		B.mutiply(vec.trans(),C);
		if(!B.eig_jacob(_eig,vec,OTHER)) return false;
		for(size_t i=1;i<=num;i++)
		{
			std::pair<Type,size_t> val_pos;
			val_pos=_eig.Lmaxa(1,i,_eig.row());
			_eig.Rswap(val_pos.second,i);
			vec.Lswap(val_pos.second,i);
		}
		vec.Rmutiply(C);
		vec.orthorized();
	}while((eig-_eig).norm()>precision);
	eig=_eig;
	return true;
}

template<class Type>
template<class U>
inline int Cmatrix<Type>::solve(const Cmatrix<U>& _rval,Cmatrix<Type>&dest)
{
	if(IsEmpty()) return -1;
	if(_rval.line()!=1) return -1;  //-1:失败;0无解;1唯一解;2无穷解
	if(row_count!=_rval.row()) return -1;
	int rtn;
	const size_t number=line_count;
	Cmatrix<size_t>index,tmp;
	combine(_rval);
	gaussed();
	Rfinda(index,0.01,GT,std::make_pair(1,row_count),std::make_pair(1,line_count));	
	index.Lfind(tmp,0,GT,std::make_pair(1,index.line()),std::make_pair(1,index.row()),OTHER);
	//tmp为秩,index为阶梯所在的列
	if(tmp(1,1)!=0)
	{
		if(index(tmp(1,1),1)==number+1) return 0;
		if(tmp(1,1)==number)
		{
			rtn=1;
			dest.resize(number,1,0);
			for(size_t i=number;i>=1;i--)
			{
				dest(i,1)=(*this)(i,number+1);
				for(size_t k=number;k>i;k--)
					dest(i,1)-=(*this)(i,k)*dest(k,1);
				dest(i,1)/=(*this)(i,i);
			}
		}
		else 
		{
			rtn=2;
			dest.resize(number,number-tmp(1,1)+1,0);
			size_t num=1;
			for(size_t i=1;i<index(1,1);i++)
			dest(i,num++)=1;
			for (size_t i=1;i<=tmp(1,1)-1;i++)
			{
				for(size_t j=index(i,1)+1;j<=index(i+1,1)-1;j++)
					dest(j,num++)=1;
			}
			for(size_t i=index(tmp(1,1),1)+1;i<=number;i++)
				dest(i,num++)=1;
			for(size_t i=1;i<=dest.line()-1;i++)
				for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),i)=0;
					for(size_t k=number;k>index(j,1);k--)
						dest(index(j,1),i)-=(*this)(j,k)*dest(k,i);
					dest(index(j,1),i)/=(*this)(j,index(j,1));
				}
			for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),dest.line())=(*this)(j,number+1);
					for(size_t k=number;k>index(j,1);k--)
						dest(index(j,1),dest.line())-=(*this)(j,k)*dest(k,dest.line());
					dest(index(j,1),dest.line())/=(*this)(j,index(j,1));
				}
		}
	}
	else
	{
		dest.resize(number,number+1,0);
		for(size_t i=1;i<=number;i++)
			dest(i,i)=1;
		rtn=2;
	}
	return rtn;
};

template<class Type>
template<class InputIterator>
inline int Cmatrix<Type>::solve(InputIterator beg,InputIterator end,Cmatrix<Type>&dest)
{
	if(IsEmpty()) return -1;
	int rtn;
	const size_t number=line_count;
	Cmatrix<size_t>index,tmp;
	push_back(beg,end);
	gaussed();
	Rfinda(index,0.01,GT,std::make_pair(1,row_count),std::make_pair(1,line_count));	
	index.Lfind(tmp,0,GT,std::make_pair(1,index.line()),std::make_pair(1,index.row()),OTHER);
	//tmp为秩,index为阶梯所在的列
	if(tmp(1,1)!=0)
	{
		if(index(tmp(1,1),1)==number+1) return 0;
		if(tmp(1,1)==number)
		{
			rtn=1;
			dest.resize(number,1,0);
			for(size_t i=number;i>=1;i--)
			{
				dest(i,1)=(*this)(i,number+1);
				for(size_t k=number;k>i;k--)
					dest(i,1)-=(*this)(i,k)*dest(k,1);
				dest(i,1)/=(*this)(i,i);
			}
		}
		else 
		{
			rtn=2;
			dest.resize(number,number-tmp(1,1)+1,0);
			size_t num=1;
			for(size_t i=1;i<index(1,1);i++)
			dest(i,num++)=1;
			for (size_t i=1;i<=tmp(1,1)-1;i++)
			{
				for(size_t j=index(i,1)+1;j<=index(i+1,1)-1;j++)
					dest(j,num++)=1;
			}
			for(size_t i=index(tmp(1,1),1)+1;i<=number;i++)
				dest(i,num++)=1;
			for(size_t i=1;i<=dest.line()-1;i++)
				for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),i)=0;
					for(size_t k=number;k>index(j,1);k--)
						dest(index(j,1),i)-=(*this)(j,k)*dest(k,i);
					dest(index(j,1),i)/=(*this)(j,index(j,1));
				}
			for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),dest.line())=(*this)(j,number+1);
					for(size_t k=number;k>index(j,1);k--)
						dest(index(j,1),dest.line())-=(*this)(j,k)*dest(k,dest.line());
					dest(index(j,1),dest.line())/=(*this)(j,index(j,1));
				}
		}
	}
	else
	{
		dest.resize(number,number+1,0);
		for(size_t i=1;i<=number;i++)
			dest(i,i)=1;
		rtn=2;
	}
	return rtn;
}

template<class Type>
template<class U>
inline int Cmatrix<Type>::solve(const U * beg,size_t num,Cmatrix<Type>&dest)
{
	if(IsEmpty()) return -1;
	int rtn;
	const size_t number=line_count;
	Cmatrix<size_t>index,tmp;
	push_back(beg,num);
	gaussed();
	Rfinda(index,0.01,GT,std::make_pair(1,row_count),std::make_pair(1,line_count));	
	index.Lfind(tmp,0,GT,std::make_pair(1,index.line()),std::make_pair(1,index.row()),OTHER);
	//tmp为秩,index为阶梯所在的列
	if(tmp(1,1)!=0)
	{
		if(index(tmp(1,1),1)==number+1) return 0;
		if(tmp(1,1)==number)
		{
			rtn=1;
			dest.resize(number,1,0);
			for(size_t i=number;i>=1;i--)
			{
				dest(i,1)=(*this)(i,number+1);
				for(size_t k=number;k>i;k--)
					dest(i,1)-=(*this)(i,k)*dest(k,1);
				dest(i,1)/=(*this)(i,i);
			}
		}
		else 
		{
			rtn=2;
			dest.resize(number,number-tmp(1,1)+1,0);
			size_t num=1;
			for(size_t i=1;i<index(1,1);i++)
			dest(i,num++)=1;
			for (size_t i=1;i<=tmp(1,1)-1;i++)
			{
				for(size_t j=index(i,1)+1;j<=index(i+1,1)-1;j++)
					dest(j,num++)=1;
			}
			for(size_t i=index(tmp(1,1),1)+1;i<=number;i++)
				dest(i,num++)=1;
			for(size_t i=1;i<=dest.line()-1;i++)
				for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),i)=0;
					for(size_t k=number;k>index(j,1);k--)
						dest(index(j,1),i)-=(*this)(j,k)*dest(k,i);
					dest(index(j,1),i)/=(*this)(j,index(j,1));
				}
			for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),dest.line())=(*this)(j,number+1);
					for(size_t k=number;k>index(j,1);k--)
						dest(index(j,1),dest.line())-=(*this)(j,k)*dest(k,dest.line());
					dest(index(j,1),dest.line())/=(*this)(j,index(j,1));
				}
		}
	}
	else
	{
		dest.resize(number,number+1,0);
		for(size_t i=1;i<=number;i++)
			dest(i,i)=1;
		rtn=2;
	}
	return rtn;
}

template<class Type>
template<class InputIterator>
inline int Cmatrix<Type>::solve(InputIterator beg,Cmatrix<Type>&dest)
{
	if(IsEmpty()) return -1;
	int rtn;
	const size_t number=line_count;
	Cmatrix<size_t>index,tmp;
	push_back(beg);
	gaussed();
	Rfinda(index,0.01,GT,std::make_pair(1,row_count),std::make_pair(1,line_count));	
	index.Lfind(tmp,0,GT,std::make_pair(1,index.line()),std::make_pair(1,index.row()),OTHER);
	//tmp为秩,index为阶梯所在的列
	if(tmp(1,1)!=0)
	{
		if(index(tmp(1,1),1)==number+1) return 0;
		if(tmp(1,1)==number)
		{
			rtn=1;
			dest.resize(number,1,0);
			for(size_t i=number;i>=1;i--)
			{
				dest(i,1)=(*this)(i,number+1);
				for(size_t k=number;k>i;k--)
					dest(i,1)-=(*this)(i,k)*dest(k,1);
				dest(i,1)/=(*this)(i,i);
			}
		}
		else 
		{
			rtn=2;
			dest.resize(number,number-tmp(1,1)+1,0);
			size_t num=1;
			for(size_t i=1;i<index(1,1);i++)
			dest(i,num++)=1;
			for (size_t i=1;i<=tmp(1,1)-1;i++)
			{
				for(size_t j=index(i,1)+1;j<=index(i+1,1)-1;j++)
					dest(j,num++)=1;
			}
			for(size_t i=index(tmp(1,1),1)+1;i<=number;i++)
				dest(i,num++)=1;
			for(size_t i=1;i<=dest.line()-1;i++)
				for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),i)=0;
					for(size_t k=number;k>index(j,1);k--)
						dest(index(j,1),i)-=(*this)(j,k)*dest(k,i);
					dest(index(j,1),i)/=(*this)(j,index(j,1));
				}
			for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),dest.line())=(*this)(j,number+1);
					for(size_t k=number;k>index(j,1);k--)
						dest(index(j,1),dest.line())-=(*this)(j,k)*dest(k,dest.line());
					dest(index(j,1),dest.line())/=(*this)(j,index(j,1));
				}
		}
	}
	else
	{
		dest.resize(number,number+1,0);
		for(size_t i=1;i<=number;i++)
			dest(i,i)=1;
		rtn=2;
	}
	return rtn;
}

template<class Type>
inline bool Cmatrix<Type>::solve_g(const Cmatrix<Type>& _rval,Cmatrix<Type>&dest)
{
	if(IsEmpty()) return false;
	if(row_count!=_rval.row()) return false;
	if(_rval.line()==0) return false;
	if(row_count!=line_count) return false;
	Cmatrix<size_t>index,tmp;
	const size_t num=line_count;
	combine(_rval);
	gaussed();
	Rfinda(index,0.01,GT,std::make_pair(1,row_count),std::make_pair(1,line_count));	
	index.Lfind(tmp,0,GT,std::make_pair(1,index.line()),std::make_pair(1,index.row()),OTHER);
	//tmp为秩,index为阶梯所在的列
	if(tmp(1,1)!=num) return false;
	dest.resize(num,_rval.line(),0);
	for(size_t j=1;j<=_rval.line();j++)
	{
		for(size_t i=num;i>=1;i--)
		{
			dest(i,j)=(*this)(i,num+j);
			for(size_t k=num;k>i;k--)
				dest(i,j)-=(*this)(i,k)*dest(k,j);
				dest(i,j)/=(*this)(i,i);
		}
	}
	return true;
};

template<class Type>
inline bool Cmatrix<Type>::solve_dolt(const Cmatrix<Type> &_rval,Cmatrix<Type>&dest)
{
	if(IsEmpty()) return false;
	if(row_count!=_rval.row()) return false;
	if(_rval.line()==0) return false;
	if(row_count!=line_count) return false;
	Cmatrix<Type> rval(_rval);
	dest.resize(row_count,rval.line());
	std::pair<Type,size_t> prim_elem;
	prim_elem=Lmaxa(1,1,row_count);
	Rswap(prim_elem.second,1);
	rval.Rswap(prim_elem.second,1);
	if(aufunc::abs((*this)(1,1))<precision) return false;
	for(size_t j=2;j<=line_count;j++)
	{
		(*this)(j,1)/=(*this)(1,1);
	}
	for(size_t i=2;i<=row_count;i++)
	{
		Cmatrix<Type> buffer(row_count+1-i,1);
		for(size_t k=1;k<=buffer.row();k++)
		{
			buffer(k,1)=(*this)(i+k-1,i);
			for(size_t j=1;j<i;j++)
				buffer(k,1)-=(*this)(i+k-1,j)*(*this)(j,i);
		}
		prim_elem=buffer.Lmaxa(1,1,buffer.row());
		Rswap(prim_elem.second+i-1,i);
		rval.Rswap(prim_elem.second+i-1,i);
		for(size_t j=1;j<i;j++)
			(*this)(i,i)-=(*this)(i,j)*(*this)(j,i);
		if(aufunc::abs((*this)(i,i))<precision) return false;
		for(size_t j=i+1;j<=line_count;j++)
		{
				//先求R的行元素
			for(size_t k=1;k<i;k++)
			{
				(*this)(i,j)-=(*this)(i,k)*(*this)(k,j);
				(*this)(j,i)-=(*this)(j,k)*(*this)(k,i);
			}
			(*this)(j,i)/=(*this)(i,i);
		}
	}
//求解方程
	for(size_t i=1;i<=rval.line();i++)
	{
		Cmatrix<Type> y(row_count,1);
		for(size_t j=1;j<=row_count;j++)
		{
			y(j,1)=rval(j,i);
			for(size_t k=1;k<j;k++)
				y(j,1)-=y(k,1)*(*this)(j,k);
		}
		for(size_t j=row_count;j>=1;j--)
		{
			dest(j,i)=y(j,1);
			for(size_t k=j+1;k<=line_count;k++)
				dest(j,i)-=dest(k,i)*(*this)(j,k);
			dest(j,i)/=(*this)(j,j);
		}
	}
	return true;
};

template<class Type>
inline bool Cmatrix<Type>::solve_chlsky(const Cmatrix<Type>& rval,Cmatrix<Type>&dest)
{
	if(!issym()) return false;
	if(row_count!=rval.row()) return false;
	if(rval.line()==0) return false;
	dest.resize(row_count,rval.line());
	for(size_t i=1;i<=row_count;i++)
	{
		for(size_t j=1;j<i;j++)
		{
			for(size_t k=1;k<j;k++)
				(*this)(i,j)-=(*this)(i,k)*(*this)(k,k)*(*this)(j,k);
			(*this)(i,j)/=(*this)(j,j);
			(*this)(i,i)-=(*this)(i,j)*(*this)(i,j)*(*this)(j,j);
		}
		if(aufunc::abs((*this)(i,i))<precision) return false;
	}
//求解方程
	for(size_t i=1;i<=rval.line();i++)
	{
		Cmatrix<Type> y(row_count,1);
		for(size_t j=1;j<=row_count;j++)
		{
			y(j,1)=rval(j,i);
			for(size_t k=1;k<j;k++)
				y(j,1)-=y(k,1)*(*this)(j,k);
		}
		for(size_t j=row_count;j>=1;j--)
		{
			dest(j,i)=y(j,1)/(*this)(j,j);
			for(size_t k=j+1;k<=line_count;k++)
				dest(j,i)-=dest(k,i)*(*this)(k,j);
		}
	}
	return true;
};

template<class Type>
inline bool Cmatrix<Type>::solve_e(const Cmatrix<Type> &_rval,Cmatrix<Type>&dest)
{
	if(issym()) return(solve_chlsky(_rval,dest));
	else return(solve_dolt(_rval,dest));
};

// 非变动形矩阵运算************************************************

//是否对称
template<class Type>
inline bool Cmatrix<Type>::issym() const
{
	using aufunc::abs;
	using std::max;
	if(IsEmpty()) return false;
	if (row_count!=line_count)
		return false;
	for (size_t i=0;i<row_count;i++)
		for (size_t j=0;j<i;j++)
			if (abs(matrix[j][i]-matrix[i][j])>precision*max(aufunc::abs(matrix[j][i]),abs(matrix[j][i])))
			return false;
	return true;
}

//是否hemit矩阵
template<class Type>
inline bool Cmatrix<Type>::ishmt() const
{
	if(IsEmpty()) return false;
	if (row_count!=line_count)
		return false;
	for (size_t i=0;i<row_count;i++)
		for (size_t j=0;j<i;j++)
			if (aufunc::abs(matrix[j][i]-aufunc::conj(matrix[i][j]))>precision*std::max(aufunc::abs(matrix[j][i]),aufunc::abs(matrix[j][i])))
			return false;
	return true;
}


//求秩
template<class Type>
inline size_t Cmatrix<Type>::rank() const
{
	if(IsEmpty()) throw bad_oper("rank出错:矩阵为空!");
	Cmatrix<Type> tmp(*this);
	std::pair<Type,size_t> vmax_pos;
	size_t pos=1; //记录搜寻主元的起始行
	for (size_t i=1;i<=tmp.line();i++)
	{
		if (pos>row_count) return pos-1;
		vmax_pos=tmp.Lmaxa(i,pos,tmp.row());
		if (aufunc::abs(vmax_pos.first)<precision) continue;
		tmp.Rswap(pos,vmax_pos.second);
		for(size_t j=pos+1;j<=tmp.row();j++)
			tmp.Radd(pos,j,j,-tmp(j,i)/tmp(pos,i));
		++pos;
	}
	return pos-1;
}


 template<class Type>
 inline Cmatrix<Type> Cmatrix<Type>::trans() const
 {
	 Cmatrix<Type> rtn(*this);
	 rtn.transed();
	 return rtn;
 };

 template<class Type>
 inline Cmatrix<Type> Cmatrix<Type>::hermit() const
 {
	 Cmatrix<Type> rtn(*this);
	 rtn.hermited();
	 return rtn;
 };

  template<class Type>
 inline Cmatrix<Type> Cmatrix<Type>::inver() const
 {
	if(IsEmpty()) throw bad_oper("inver出错:矩阵为空!");
	if(line_count!=row_count) return *this;
	Cmatrix<Type> tmp(row_count,line_count*2);
	tmp.united(row_count);
	tmp.combine(*this);
	std::pair<Type,size_t> vmax_pos;
	size_t pos=1; //记录搜寻主元的起始行
	for (size_t i=line_count+1;i<=line_count*2;i++)
	{
		vmax_pos=tmp.Lmaxa(i,pos,tmp.row());
		if (aufunc::abs(vmax_pos.first)<precision) return *this;
		tmp.Rswap(pos,vmax_pos.second);
		tmp.Rscale(pos,tmp(pos,i),false);
		for(size_t j=1;j<=tmp.row();j++)
		{
			if(j==pos) continue;
			tmp.Radd(pos,j,j,-tmp(j,i));
		}
		++pos;
	}
	tmp.Lresize(line_count);
	return tmp;
 }

template<class Type>
inline Cmatrix<Type>  Cmatrix<Type>::rever(Method flag=LINE) const
{
	Cmatrix<Type> tmp(*this);
	tmp.revered(flag);
	return tmp;
};

template<class Type>
inline Cmatrix<Type>  Cmatrix<Type>::std() const
{
	Cmatrix<Type> tmp(*this);
	tmp.stded();
	return tmp;
}

template<class Type>
inline Cmatrix<Type>  Cmatrix<Type>::gauss() const
{
	Cmatrix<Type> tmp(*this);
	tmp.gaussed();
	return tmp;
}

template<class Type>
inline void Cmatrix<Type>::fft(Cmatrix<std::complex<Type> >&dest,Method flag=LINE) const
{
	if(IsEmpty()) throw bad_oper("fft出错:矩阵为空!");
	dest.resize(row_count,line_count,std::complex<Type>(0,0));
	const long double PI=3.1415926;
	if(flag==LINE)
	{
		const size_t N=row_count;
		for(size_t i=1;i<=line_count;i++)
		{
			std::vector<std::complex<Type> > buffer(N);
			for(size_t j=0;j<N;++j)
				buffer[j]=std::exp(std::complex<Type>(0,-2*PI*j/N));
			for(size_t j=1;j<=N;j++)
				for(size_t k=0;k<N;k++)
					dest(j,i)=dest(j,i)+buffer[((j-1)*k)%N]*(*this)(k+1,i);
		}		
	}
	else
	{
		const size_t N=line_count;
		for(size_t i=1;i<=row_count;i++)
		{
			std::vector<std::complex<Type> > buffer(N);
			for(size_t j=0;j<N;++j)
				buffer[j]=std::exp(std::complex<Type>(0,-2*PI*j/N));
			for(size_t j=1;j<=N;j++)
				for(size_t k=0;k<N;k++)
					dest(i,j)=dest(i,j)+buffer[((j-1)*k)%N]*(*this)(i,k+1);
		}
	}
	
}

template<class Type>
inline void Cmatrix<Type>::fft(Cmatrix<Type>&dest,Method flag=LINE) const
{
	if(IsEmpty()) throw bad_oper("fft出错:矩阵为空!");
	dest.resize(row_count,line_count,Type(0,0));
	const long double PI=3.1415926;
	if(flag==LINE)
	{
		const size_t N=row_count;
		for(size_t i=1;i<=line_count;i++)
		{
			std::vector<Type> buffer(N);
			for(size_t j=0;j<N;++j)
				buffer[j]=std::exp(Type(0,-2*PI*j/N));
			for(size_t j=1;j<=N;j++)
				for(size_t k=0;k<N;k++)
					dest(j,i)=dest(j,i)+buffer[((j-1)*k)%N]*(*this)(k+1,i);
		}
		
	}
	else
	{
		const size_t N=line_count;
		for(size_t i=1;i<=row_count;i++)
		{
			std::vector<Type> buffer(N);
			for(size_t j=0;j<N;++j)
				buffer[j]=std::exp(Type(0,-2*PI*j/N));
			for(size_t j=1;j<=N;j++)
				for(size_t k=0;k<N;k++)
					dest(i,j)=dest(i,j)+buffer[((j-1)*k)%N]*(*this)(i,k+1);
		}
		
	}
}

template<class Type>
inline void Cmatrix<Type>::fftn(Cmatrix<std::complex<Type> >&dest,Method flag=LINE) const
{
	if(IsEmpty()) throw bad_oper("fftn出错:矩阵为空!");
	dest.resize(row_count,line_count,std::complex<Type>(0,0));
	const long double PI=3.1415926;
	if(flag==LINE)
	{
		const size_t N=row_count;
		for(size_t i=1;i<=line_count;i++)
		{
			std::vector<std::complex<Type> > buffer(N);
			for(size_t j=0;j<N;++j)
				buffer[j]=std::exp(std::complex<Type>(0,2*PI*j/N));
			for(size_t j=1;j<=N;j++)
				for(size_t k=0;k<N;k++)
					dest(j,i)=dest(j,i)+buffer[((j-1)*k)%N]*(*this)(k+1,i)/Type(N);
		}		
	}
	else
	{
		const size_t N=line_count;
		for(size_t i=1;i<=row_count;i++)
		{
			std::vector<std::complex<Type> > buffer(N);
			for(size_t j=0;j<N;++j)
				buffer[j]=std::exp(std::complex<Type>(0,2*PI*j/N));
			for(size_t j=1;j<=N;j++)
				for(size_t k=0;k<N;k++)
					dest(i,j)=dest(i,j)+buffer[((j-1)*k)%N]*(*this)(i,k+1)/Type(N);
		}
	}
	
}

template<class Type>
inline void Cmatrix<Type>::fftn(Cmatrix<Type>&dest,Method flag=LINE) const
{
	if(IsEmpty()) throw bad_oper("fftn出错:矩阵为空!");
	dest.resize(row_count,line_count,Type(0,0));
	const long double PI=3.1415926;
	if(flag==LINE)
	{
		const size_t N=row_count;
		for(size_t i=1;i<=line_count;i++)
		{
			std::vector<Type> buffer(N);
			for(size_t j=0;j<N;++j)
				buffer[j]=std::exp(Type(0,2*PI*j/N));
			for(size_t j=1;j<=N;j++)
				for(size_t k=0;k<N;k++)
					dest(j,i)=dest(j,i)+buffer[((j-1)*k)%N]*(*this)(k+1,i)/Type(N,0);
		}
		
	}
	else
	{
		const size_t N=line_count;
		for(size_t i=1;i<=row_count;i++)
		{
			std::vector<Type> buffer(N);
			for(size_t j=0;j<N;++j)
				buffer[j]=std::exp(Type(0,-2*PI*j/N));
			for(size_t j=1;j<=N;j++)
				for(size_t k=0;k<N;k++)
					dest(i,j)=dest(i,j)+buffer[((j-1)*k)%N]*(*this)(i,k+1)/Type(N,0);
		}
		
	}
}

template<class Type>
inline Type Cmatrix<Type>::mode() const
{
	if(IsEmpty()) throw bad_oper("mode出错:矩阵为空!");
	if(row_count!=line_count) throw  bad_oper("矩阵不为方阵!");
	Cmatrix<Type> tmp(*this);
	std::pair<Type,size_t> vmax_pos;
	size_t pos=1; //记录搜寻主元的起始行
	Type rtn=Type(1);
	for (size_t i=1;i<=tmp.line();i++)
	{
		vmax_pos=tmp.Lmaxa(i,pos,tmp.row());
		if (aufunc::abs(vmax_pos.first)<precision) return Type(0);
		tmp.Rswap(pos,vmax_pos.second);
		for(size_t j=pos+1;j<=tmp.row();j++)
			tmp.Radd(pos,j,j,-tmp(j,i)/tmp(pos,i));
		++pos;
		rtn=rtn*tmp(i,i);
	}
	return rtn;
}

template<class Type>
inline int Cmatrix<Type>::quadform() const
{
	if(line_count!=row_count||!issym()) return 0;
	Type tmp=getzone(1,1,1,1).mode();
	int rtn=0;
	for (size_t i=1;i<=line_count;i++)
	{
		Type tmp=getzone(1,1,i,i).mode();
		if(i%2==1)
		{
			if(tmp>0)
			{
				if(rtn<0) return 0;
				else if(rtn==0) rtn=1;
			}
			else if(tmp<0)
			{
				if(rtn>0) return 0;
				else if(rtn==0) rtn=-1;
			}
			else
			{
				if(rtn==1) rtn=2;
				else if(rtn==-1) rtn=-2;
			}
		}
		else
		{
			if(tmp<0) return 0;
			else if(tmp==0)
			{
				if(rtn=-1) rtn=-2;
				else if(rtn=1) rtn=2;
			}
		}
	}
	return rtn;
}

template<class Type>
inline bool Cmatrix<Type>::dolt(Cmatrix<Type>& L,Cmatrix<Type>&R,Method flag=OTHER) const
{
	if(IsEmpty()) return false;
	if(row_count!=line_count) return false;
	L.resize(row_count,line_count);
	R.resize(row_count,line_count);
	L=0;
	R=0;
	if(flag==LZY)
	{
		Cmatrix<Type> tmp(*this);
		L(1,1)=Type(1);
		std::pair<Type,size_t> prim_elem;
		prim_elem=tmp.Lmaxa(1,1,tmp.row());
		tmp.Rswap(prim_elem.second,1);
		R(1,1)=tmp(1,1);
		if(aufunc::abs(R(1,1)<precision)) return false;
		for(size_t j=2;j<=line_count;j++)
		{
			R(1,j)=tmp(1,j)/L(1,1);
			L(j,1)=tmp(j,1)/R(1,1);
		}
		for(size_t i=2;i<=row_count;i++)
		{
			Cmatrix<Type> buffer(row_count+1-i);
			for(size_t k=1;k<=buffer.row();k++)
			{
				buffer(k,1)=tmp(i+k-1,i);
				for(size_t j=1;j<i;j++)
					buffer(k,1)-=L(i+k-1,j)*R(j,i);
			}
			prim_elem=buffer.Lmaxa(1,1,buffer.row());
			tmp.Rswap(prim_elem.second+i-1,i);
			L.Rswap(prim_elem.second+i-1,i);
			L(i,i)=Type(1);
			R(i,i)=tmp(i,i);
			for(size_t j=1;j<i;j++)
				R(i,i)-=L(i,j)*R(j,i);
			if(aufunc::abs(R(i,i))<precision) return false;
			for(size_t j=i+1;j<=line_count;j++)
			{
				//先求R的行元素
				R(i,j)=tmp(i,j);
				L(j,i)=tmp(j,i);
				for(size_t k=1;k<i;k++)
				{
					R(i,j)-=L(i,k)*R(k,j);
					L(j,i)-=L(j,k)*R(k,i);
				}
				L(j,i)/=R(i,i);
			}
		}
	}
	else
	{
		L(1,1)=Type(1);
		R(1,1)=(*this)(1,1);
		if(aufunc::abs(R(1,1))<precision) return false;
		for(size_t j=2;j<=line_count;j++)
		{
			R(1,j)=(*this)(1,j)/L(1,1);
			L(j,1)=(*this)(j,1)/R(1,1);
		}
		for(size_t i=2;i<=row_count;i++)
		{
			L(i,i)=Type(1);
			R(i,i)=(*this)(i,i);
			for(size_t j=1;j<i;j++)
				R(i,i)-=L(i,j)*R(j,i);
			if(aufunc::abs(R(i,i))<precision) return false;
			for(size_t j=i+1;j<=line_count;j++)
			{
				//先求R的行元素
				R(i,j)=(*this)(i,j);
				L(j,i)=(*this)(j,i);
				for(size_t k=1;k<i;k++)
				{
					R(i,j)-=L(i,k)*R(k,j);
					L(j,i)-=L(j,k)*R(k,i);
				}
				L(j,i)/=R(i,i);
			}
		}
	}
	return true;
}

template<class Type>
inline bool Cmatrix<Type>::crout(Cmatrix<Type>& L,Cmatrix<Type>&R,Method flag=OTHER) const
{
	if(IsEmpty()) return false;
	if(row_count!=line_count) return false;
	L.resize(row_count,line_count);
	R.resize(row_count,line_count);
	L=0;
	R=0;
	if(flag==LZY)
	{
		Cmatrix<Type> tmp(*this);
		R(1,1)=Type(1);
		std::pair<Type,size_t> prim_elem;
		prim_elem=tmp.Lmaxa(1,1,tmp.row());
		tmp.Rswap(prim_elem.second,1);
		L(1,1)=tmp(1,1);
		if(aufunc::abs(L(1,1))<precision) return false;
		for(size_t j=2;j<=line_count;j++)
		{
			R(1,j)=tmp(1,j)/L(1,1);
			L(j,1)=tmp(j,1)/R(1,1);
		}
		for(size_t i=2;i<=row_count;i++)
		{
			Cmatrix<Type> buffer(row_count+1-i);
			for(size_t k=1;k<=buffer.row();k++)
			{
				buffer(k,1)=tmp(i+k-1,i);
				for(size_t j=1;j<i;j++)
					buffer(k,1)-=L(i+k-1,j)*R(j,i);
			}
			prim_elem=buffer.Lmaxa(1,1,buffer.row());
			tmp.Rswap(prim_elem.second+i-1,i);
			L.Rswap(prim_elem.second+i-1,i);
			R(i,i)=Type(1);
			L(i,i)=tmp(i,i);
			for(size_t j=1;j<i;j++)
				L(i,i)-=L(i,j)*R(j,i);
			if(aufunc::abs(L(i,i))<precision) return false;
			for(size_t j=i+1;j<=line_count;j++)
			{
				R(i,j)=tmp(i,j);
				L(j,i)=tmp(j,i);
				for(size_t k=1;k<i;k++)
				{
					R(i,j)-=L(i,k)*R(k,j);
					L(j,i)-=L(j,k)*R(k,i);
				}
				R(i,j)/=L(i,i);
			}
		}
	}
	else
	{
		R(1,1)=Type(1);
		L(1,1)=(*this)(1,1);
		if(aufunc::abs(L(1,1))<precision) return false;
		for(size_t j=2;j<=line_count;j++)
		{
			R(1,j)=(*this)(1,j)/L(1,1);
			L(j,1)=(*this)(j,1)/R(1,1);
		}
		for(size_t i=2;i<=row_count;i++)
		{
			R(i,i)=Type(1);
			L(i,i)=(*this)(i,i);
			for(size_t j=1;j<i;j++)
				L(i,i)-=L(i,j)*R(j,i);
			if(aufunc::abs(L(i,i))<precision) return false;
			for(size_t j=i+1;j<=line_count;j++)
			{
				R(i,j)=(*this)(i,j);
				L(j,i)=(*this)(j,i);
				for(size_t k=1;k<i;k++)
				{
					R(i,j)-=L(i,k)*R(k,j);
					L(j,i)-=L(j,k)*R(k,i);
				}
				R(i,j)/=L(i,i);
			}
		}
	}
	return true;
}

template<class Type>
inline bool Cmatrix<Type>::chlsky(Cmatrix<Type>& L,Cmatrix<Type>&D) const
{
	if(!issym()) return false;
	D.resize(row_count,line_count);
	L.united(row_count);
	D=0;
	D(1,1)=(*this)(1,1);
	if(aufunc::abs(D(1,1))<precision)return false;
	for(size_t i=2;i<row_count;i++)
	{
		D(i,i)=(*this)(i,i);
		for(size_t j=1;j<i;j++)
		{
			L(i,j)=(*this)(i,j);
			for(size_t k=1;k<j;k++)
				L(i,j)-=L(i,k)*D(k,k)*L(j,k);
			L(i,j)/=D(j,j);
			D(i,i)-=L(i,j)*L(i,j)*D(j,j);
		}
		if(aufunc::abs(D(i,i))<precision) return false;
	}
	D(row_count,row_count)=(*this)(row_count,row_count);
	for(size_t j=1;j<row_count;j++)
	{
		L(row_count,j)=(*this)(row_count,j);
		for(size_t k=1;k<j;k++)
			L(row_count,j)-=L(row_count,k)*D(k,k)*L(j,k);
		L(row_count,j)/=D(j,j);
		D(row_count,row_count)-=L(row_count,j)*L(row_count,j)*D(j,j);	
	}
	return true;
}

template<class Type>
inline Cmatrix<Type >Cmatrix<Type>::zero() const
{
	Cmatrix<Type> rst(*this);
	rst.zeroed();
	return rst;
}

template<class Type>
inline long double Cmatrix<Type>::radius(Type orig=Type(),const size_t num=100) const
{
	if(IsEmpty()) throw bad_oper("radius出错:矩阵为空!");
	if(row_count!=line_count) throw bad_oper("radius出错:矩阵不为方阵!");
	Type rtn=0;
	Cmatrix<Type> v(row_count,1,1),u(row_count,1,1);
	rtn=0;
	size_t i=0;
	size_t j=0;
	if(orig==Type())
	{
		do
		{
			u.Rmutiply((*this));
			rtn=u.Lmaxa(1,1,u.row()).first;
			u=u/rtn;
			v-=u;
			for(i=1;i<=row_count;i++)
				if(aufunc::abs(v(i,1))>precision) break;
			v=u;
			j++;
		}while(i<=row_count&&j<=num);
	}
	else
	{
		Cmatrix tmp(*this);
		for(size_t k=1;k<=row_count;k++)
			tmp(k,k)-=orig;
		do
		{
			u.Rmutiply(tmp);
			rtn=u.Lmaxa(1,1,u.row()).first;
			u=u/rtn;
			v-=u;
			for(i=1;i<=row_count;i++)
				if(aufunc::abs(v(i,1))>precision) break;
			v=u;
			j++;
		}while(i<=row_count&&j<=num);
		rtn+=orig;
	}
	if(i>row_count) return long double(aufunc::abs(rtn));
	else throw bad_oper("radius出错:求解失败!");
	
}

template<class Type>
inline bool Cmatrix<Type>::eig_jacob (Cmatrix<Type>& eig,Cmatrix<Type>&vec,Method flag=INC) const
{
	if(!issym()) return false;
	Cmatrix<Type> tmp(*this);
	eig.resize(row_count,1);
	vec.united(row_count);
	long double value=0;
	const long double n=(row_count*line_count-row_count)/2.0;
	value=0;
	for(size_t i=1;i<=row_count;i++)
		for(size_t j=1;j<i;j++)
			value+=aufunc::abs(tmp(i,j)*tmp(i,j));
	do
	{
		//计算非对角元素的平方和
		value/=n;
		size_t lpos=0,rpos=0;
		for(rpos=1;rpos<=row_count;rpos++)
		{
			for(lpos=rpos+1;lpos<=line_count;lpos++)
			{
				if(aufunc::abs(tmp(rpos,lpos))>value)
				{
					long double y=aufunc::abs(tmp(rpos,rpos)-tmp(lpos,lpos));
					long double x=2.0*tmp(rpos,lpos);
					if(tmp(rpos,rpos)<tmp(lpos,lpos)) x*=(-1.0);
					long double cos2a=y/std::sqrt(x*x+y*y);
					long double cosa=std::sqrt((1+cos2a)/2);
					long double sin2a=x/std::sqrt(x*x+y*y);
					long double sina=sin2a/(2*cosa);
					const long double pq=tmp(rpos,lpos),pp=tmp(rpos,rpos),qq=tmp(lpos,lpos);
					for(size_t j=1;j<=line_count;j++)
					{
						long double qj=tmp(lpos,j),pj=tmp(rpos,j);
						tmp(j,rpos)=tmp(rpos,j)=pj*cosa+qj*sina;
						tmp(lpos,j)=tmp(j,lpos)=(0-pj)*sina+qj*cosa;
					}
					tmp(rpos,rpos)=pp*cosa*cosa+2*pq*sina*cosa+qq*sina*sina;
					tmp(lpos,lpos)=pp*sina*sina-2*pq*sina*cosa+qq*cosa*cosa;
					tmp(rpos,lpos)=tmp(lpos,rpos)=0.5*(qq-pp)*sin2a+pq*cos2a;
					Cmatrix<Type> R(UNIT,row_count);
					R(rpos,rpos)=R(lpos,lpos)=cosa;
					R(rpos,lpos)=-sina;R(lpos,rpos)=sina;
					vec*=R;
				}
			}
		}
	}while(value>precision);
	for(size_t i=1;i<=row_count;i++)
		eig(i,1)=tmp(i,i);
	if(flag==INC)
	{
		for(size_t i=1;i<=row_count;i++)
		{
			std::pair<Type,size_t> val_pos;
			val_pos=eig.Lmin(1,i,eig.row());
			eig.Rswap(val_pos.second,i);
			vec.Lswap(val_pos.second,i);
		}
	}
	else if(flag==DEC)
	{
		for(size_t i=1;i<=row_count;i++)
		{
			std::pair<Type,size_t> val_pos;
			val_pos=eig.Lmax(1,i,eig.row());
			eig.Rswap(val_pos.second,i);
			vec.Lswap(val_pos.second,i);
		}
	}
	return true;
}

template<class Type>
inline bool Cmatrix<Type>::eig_subspace(Cmatrix<Type>& eig,Cmatrix<Type>& vec,const size_t num=1) const
{
	if(!issym()) return false;
	if(num>line_count) return false;
	eig.resize(num,1);
	vec.resize(row_count,num);
	Cmatrix<Type> C(row_count,num),B(num,num);
	for(size_t i=1;i<=num;i++)
		vec(i,i)=1;
	vec.orthorized();
	Cmatrix<Type> _eig(num);
	do
	{
		eig=_eig;
		C=vec;
		C.Rmutiply(*this);
		B.mutiply(vec.trans(),C);
		if(!B.eig_jacob(_eig,vec,OTHER)) return false;
		for(size_t i=1;i<=num;i++)
		{
			std::pair<Type,size_t> val_pos;
			val_pos=_eig.Lmaxa(1,i,_eig.row());
			_eig.Rswap(val_pos.second,i);
			vec.Lswap(val_pos.second,i);
		}
		vec.Rmutiply(C);
		vec.orthorized();
	}while((eig-_eig).norm()>precision);
	eig=_eig;
	return true;
}

template<class Type>
template<class U>
inline int Cmatrix<Type>::solve(const Cmatrix<U>& _rval,Cmatrix<Type>&dest) const
{
	if(IsEmpty()) return -1;;
	if(_rval.line()!=1) return -1;  //-1:失败;0无解;1唯一解;2无穷解
	if(row_count!=_rval.row()) return -1;
	int rtn;
	Cmatrix<Type> orig(*this);
	Cmatrix<size_t>index,tmp;
	orig.combine(_rval);
	orig.gaussed();
	orig.Rfinda(index,0.01,GT,std::make_pair(1,orig.row()),std::make_pair(1,orig.line()));	
	index.Lfind(tmp,0,GT,std::make_pair(1,index.line()),std::make_pair(1,index.row()),OTHER);
	//tmp为秩,index为阶梯所在的列
	if(tmp(1,1)!=0)
	{
		if(index(tmp(1,1),1)==line_count+1) return 0;
		if(tmp(1,1)==line_count)
		{
			rtn=1;
			dest.resize(line_count,1,0);
			for(size_t i=line_count;i>=1;i--)
			{
				dest(i,1)=orig(i,line_count+1);
				for(size_t k=line_count;k>i;k--)
					dest(i,1)-=orig(i,k)*dest(k,1);
				dest(i,1)/=orig(i,i);
			}
		}
		else 
		{
			rtn=2;
			dest.resize(line_count,line_count-tmp(1,1)+1,0);
			size_t num=1;
			for(size_t i=1;i<index(1,1);i++)
			dest(i,num++)=1;
			for (size_t i=1;i<=tmp(1,1)-1;i++)
			{
				for(size_t j=index(i,1)+1;j<=index(i+1,1)-1;j++)
					dest(j,num++)=1;
			}
			for(size_t i=index(tmp(1,1),1)+1;i<=line_count;i++)
				dest(i,num++)=1;
			for(size_t i=1;i<=dest.line()-1;i++)
				for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),i)=0;
					for(size_t k=line_count;k>index(j,1);k--)
						dest(index(j,1),i)-=orig(j,k)*dest(k,i);
					dest(index(j,1),i)/=orig(j,index(j,1));
				}
			for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),dest.line())=orig(j,line_count+1);
					for(size_t k=line_count;k>index(j,1);k--)
						dest(index(j,1),dest.line())-=orig(j,k)*dest(k,dest.line());
					dest(index(j,1),dest.line())/=orig(j,index(j,1));
				}
		}
	}
	else
	{
		dest.resize(line_count,line_count+1,0);
		for(size_t i=1;i<=line_count;i++)
			dest(i,i)=1;
		rtn=2;
	}
	return rtn;
};

template<class Type>
template<class InputIterator>
inline int Cmatrix<Type>::solve(InputIterator beg,InputIterator end,Cmatrix<Type>&dest) const
{
	if(IsEmpty()) return -1;
	int rtn;
	Cmatrix<Type> orig(*this);
	Cmatrix<size_t>index,tmp;
	orig.push_back(beg,end);
	orig.gaussed();
	orig.Rfinda(index,0.01,GT,std::make_pair(1,orig.row()),std::make_pair(1,orig.line()));	
	index.Lfind(tmp,0,GT,std::make_pair(1,index.line()),std::make_pair(1,index.row()),OTHER);
	//tmp为秩,index为阶梯所在的列
	if(tmp(1,1)!=0)
	{
		if(index(tmp(1,1),1)==line_count+1) return 0;
		if(tmp(1,1)==line_count)
		{
			rtn=1;
			dest.resize(line_count,1,0);
			for(size_t i=line_count;i>=1;i--)
			{
				dest(i,1)=orig(i,line_count+1);
				for(size_t k=line_count;k>i;k--)
					dest(i,1)-=orig(i,k)*dest(k,1);
				dest(i,1)/=orig(i,i);
			}
		}
		else 
		{
			rtn=2;
			dest.resize(line_count,line_count-tmp(1,1)+1,0);
			size_t num=1;
			for(size_t i=1;i<index(1,1);i++)
			dest(i,num++)=1;
			for (size_t i=1;i<=tmp(1,1)-1;i++)
			{
				for(size_t j=index(i,1)+1;j<=index(i+1,1)-1;j++)
					dest(j,num++)=1;
			}
			for(size_t i=index(tmp(1,1),1)+1;i<=line_count;i++)
				dest(i,num++)=1;
			for(size_t i=1;i<=dest.line()-1;i++)
				for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),i)=0;
					for(size_t k=line_count;k>index(j,1);k--)
						dest(index(j,1),i)-=orig(j,k)*dest(k,i);
					dest(index(j,1),i)/=orig(j,index(j,1));
				}
			for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),dest.line())=orig(j,line_count+1);
					for(size_t k=line_count;k>index(j,1);k--)
						dest(index(j,1),dest.line())-=orig(j,k)*dest(k,dest.line());
					dest(index(j,1),dest.line())/=orig(j,index(j,1));
				}
		}
	}
	else
	{
		dest.resize(line_count,line_count+1,0);
		for(size_t i=1;i<=line_count;i++)
			dest(i,i)=1;
		rtn=2;
	}
	return rtn;
}

template<class Type>
template<class U>
inline int Cmatrix<Type>::solve(const U * beg,size_t num,Cmatrix<Type>&dest) const
{
	if(IsEmpty()) return -1;
	int rtn;
	Cmatrix<Type> orig(*this);
	Cmatrix<size_t>index,tmp;
	orig.push_back(beg,num);
	orig.gaussed();
	orig.Rfinda(index,0.01,GT,std::make_pair(1,orig.row()),std::make_pair(1,orig.line()));	
	index.Lfind(tmp,0,GT,std::make_pair(1,index.line()),std::make_pair(1,index.row()),OTHER);
	//tmp为秩,index为阶梯所在的列
	if(tmp(1,1)!=0)
	{
		if(index(tmp(1,1),1)==line_count+1) return 0;
		if(tmp(1,1)==line_count)
		{
			rtn=1;
			dest.resize(line_count,1,0);
			for(size_t i=line_count;i>=1;i--)
			{
				dest(i,1)=orig(i,line_count+1);
				for(size_t k=line_count;k>i;k--)
					dest(i,1)-=orig(i,k)*dest(k,1);
				dest(i,1)/=orig(i,i);
			}
		}
		else 
		{
			rtn=2;
			dest.resize(line_count,line_count-tmp(1,1)+1,0);
			size_t num=1;
			for(size_t i=1;i<index(1,1);i++)
			dest(i,num++)=1;
			for (size_t i=1;i<=tmp(1,1)-1;i++)
			{
				for(size_t j=index(i,1)+1;j<=index(i+1,1)-1;j++)
					dest(j,num++)=1;
			}
			for(size_t i=index(tmp(1,1),1)+1;i<=line_count;i++)
				dest(i,num++)=1;
			for(size_t i=1;i<=dest.line()-1;i++)
				for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),i)=0;
					for(size_t k=line_count;k>index(j,1);k--)
						dest(index(j,1),i)-=orig(j,k)*dest(k,i);
					dest(index(j,1),i)/=orig(j,index(j,1));
				}
			for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),dest.line())=orig(j,line_count+1);
					for(size_t k=line_count;k>index(j,1);k--)
						dest(index(j,1),dest.line())-=orig(j,k)*dest(k,dest.line());
					dest(index(j,1),dest.line())/=orig(j,index(j,1));
				}
		}
	}
	else
	{
		dest.resize(line_count,line_count+1,0);
		for(size_t i=1;i<=line_count;i++)
			dest(i,i)=1;
		rtn=2;
	}
	return rtn;
}

template<class Type>
template<class InputIterator>
inline int Cmatrix<Type>::solve(InputIterator beg,Cmatrix<Type>&dest) const
{
	if(IsEmpty()) return -1;
	int rtn;
	Cmatrix<Type> orig(*this);
	Cmatrix<size_t>index,tmp;
	orig.push_back(beg);
	orig.gaussed();
	orig.Rfinda(index,0.01,GT,std::make_pair(1,orig.row()),std::make_pair(1,orig.line()));	
	index.Lfind(tmp,0,GT,std::make_pair(1,index.line()),std::make_pair(1,index.row()),OTHER);
	//tmp为秩,index为阶梯所在的列
	if(tmp(1,1)!=0)
	{
		if(index(tmp(1,1),1)==line_count+1) return 0;
		if(tmp(1,1)==line_count)
		{
			rtn=1;
			dest.resize(line_count,1,0);
			for(size_t i=line_count;i>=1;i--)
			{
				dest(i,1)=orig(i,line_count+1);
				for(size_t k=line_count;k>i;k--)
					dest(i,1)-=orig(i,k)*dest(k,1);
				dest(i,1)/=orig(i,i);
			}
		}
		else 
		{
			rtn=2;
			dest.resize(line_count,line_count-tmp(1,1)+1,0);
			size_t num=1;
			for(size_t i=1;i<index(1,1);i++)
			dest(i,num++)=1;
			for (size_t i=1;i<=tmp(1,1)-1;i++)
			{
				for(size_t j=index(i,1)+1;j<=index(i+1,1)-1;j++)
					dest(j,num++)=1;
			}
			for(size_t i=index(tmp(1,1),1)+1;i<=line_count;i++)
				dest(i,num++)=1;
			for(size_t i=1;i<=dest.line()-1;i++)
				for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),i)=0;
					for(size_t k=line_count;k>index(j,1);k--)
						dest(index(j,1),i)-=orig(j,k)*dest(k,i);
					dest(index(j,1),i)/=orig(j,index(j,1));
				}
			for(size_t j=tmp(1,1);j>=1;j--)
				{
					dest(index(j,1),dest.line())=orig(j,line_count+1);
					for(size_t k=line_count;k>index(j,1);k--)
						dest(index(j,1),dest.line())-=orig(j,k)*dest(k,dest.line());
					dest(index(j,1),dest.line())/=orig(j,index(j,1));
				}
		}
	}
	else
	{
		dest.resize(line_count,line_count+1,0);
		for(size_t i=1;i<=line_count;i++)
			dest(i,i)=1;
		rtn=2;
	}
	return rtn;
}

template<class Type>
inline bool Cmatrix<Type>::solve_g(const Cmatrix<Type>& _rval,Cmatrix<Type>&dest) const
{
	if(IsEmpty()) return false;
	if(row_count!=_rval.row()) return false;
	if(_rval.line()==0) return false;
	if(row_count!=line_count) return false;
	Cmatrix<Type> orig(*this);
	Cmatrix<size_t>index,tmp;
	orig.combine(_rval);
	orig.gaussed();
	orig.Rfinda(index,0.01,GT,std::make_pair(1,orig.row()),std::make_pair(1,orig.line()));	
	index.Lfind(tmp,0,GT,std::make_pair(1,index.line()),std::make_pair(1,index.row()),OTHER);
	//tmp为秩,index为阶梯所在的列
	if(tmp(1,1)!=line_count) return false;
	dest.resize(line_count,_rval.line(),0);
	for(size_t j=1;j<=_rval.line();j++)
	{
		for(size_t i=line_count;i>=1;i--)
		{
			dest(i,j)=orig(i,line_count+j);
			for(size_t k=line_count;k>i;k--)
				dest(i,j)-=orig(i,k)*dest(k,j);
				dest(i,j)/=orig(i,i);
		}
	}
	return true;
};

template<class Type>
inline bool Cmatrix<Type>::solve_dolt(const Cmatrix<Type> &_rval,Cmatrix<Type>&dest) const
{
	Cmatrix<Type> tmp(*this);
	return tmp.solve_dolt(_rval,dest);
};

template<class Type>
inline bool Cmatrix<Type>::solve_chlsky(const Cmatrix<Type>& rval,Cmatrix<Type>&dest) const
{
	if(!issym()) return false;
	if(row_count!=rval.row()) return false;
	if(rval.line()==0) return false;
	Cmatrix<Type> L(UNIT,row_count);
	dest.resize(row_count,rval.line());
	for(size_t i=1;i<=row_count;i++)
	{
		L(i,i)=(*this)(i,i);
		for(size_t j=1;j<i;j++)
		{
			L(i,j)=(*this)(i,j);
			for(size_t k=1;k<j;k++)
				L(i,j)-=L(i,k)*L(k,k)*L(j,k);
			L(i,j)/=L(j,j);
			L(i,i)-=L(i,j)*L(i,j)*L(j,j);
		}
		if(aufunc::abs(L(i,i))<precision) return false;
	}
//求解方程
	for(size_t i=1;i<=rval.line();i++)
	{
		Cmatrix<Type> y(row_count,1);
		for(size_t j=1;j<=row_count;j++)
		{
			y(j,1)=rval(j,i);
			for(size_t k=1;k<j;k++)
				y(j,1)-=y(k,1)*L(j,k);
		}
		for(size_t j=row_count;j>=1;j--)
		{
			dest(j,i)=y(j,1)/L(j,j);
			for(size_t k=j+1;k<=line_count;k++)
				dest(j,i)-=dest(k,i)*L(k,j);
		}
	}
	return true;
};

template<class Type>
inline bool Cmatrix<Type>::solve_e(const Cmatrix<Type> &_rval,Cmatrix<Type>&dest) const
{
	if(issym()) return(solve_chlsky(_rval,dest));
	else return(solve_dolt(_rval,dest));
};

template<class Type>
inline long double Cmatrix<Type>::norm(Method flag) const
{
	if(IsEmpty()) throw bad_oper("norm出错:矩阵为空!");
	long double rtn=0;
	if(flag==M1)
	{
		for(size_t i=1;i<=row_count;i++)
			for(size_t j=1;j<=line_count;j++)
				rtn+=aufunc::abs((*this)(i,j));
	}
	else if(flag==F)
	{
		for(size_t i=1;i<=row_count;i++)
			for(size_t j=1;j<=line_count;j++)
				rtn+=aufunc::abs((*this)(i,j)*(*this)(i,j));
		rtn=std::pow(rtn,long double(0.5));
	}
	else if(flag==G)
	{
		rtn=aufunc::abs((*this)(1,1));
		for(size_t i=1;i<=row_count;i++)
			for(size_t j=1;j<=line_count;j++)
				if(rtn<aufunc::abs((*this)(i,j))) rtn=aufunc::abs((*this)(i,j));
		rtn*=std::max(row_count,line_count);
	}
	else
	{
		rtn=aufunc::abs((*this)(1,1));
		for(size_t i=1;i<=row_count;i++)
			for(size_t j=1;j<=line_count;j++)
				if(rtn<aufunc::abs((*this)(i,j))) rtn=aufunc::abs((*this)(i,j));
		rtn*=std::pow(long double(row_count*line_count),long double(0.5));
	}
	return rtn;
}

template<class Type>
inline long double Cmatrix<Type>::norm(int flag=1) const
{
	if(IsEmpty()) throw bad_oper("norm出错:矩阵为空!");
	long double rtn=0;
	if(flag==1)
	{
		for(size_t j=1;j<=line_count;j++)
		{
			long double tmp=0;
			for(size_t i=1;i<=row_count;i++)
				tmp+=aufunc::abs((*this)(i,j));
			if(rtn<tmp) rtn=tmp;
		}
	}
	else if(flag==2)
	{
		Cmatrix<Type> tmp(*this);
		tmp.hermited();
		tmp*=(*this);
		rtn=std::pow(long double(aufunc::abs(tmp.radius())),long double(0.5));
	}
	else
	{
		for(size_t i=1;i<=row_count;i++)
		{
			long double tmp=0;
			for(size_t j=1;j<=line_count;j++)
				tmp+=aufunc::abs((*this)(i,j));
			if(rtn<tmp) rtn=tmp;
		}
	}
	return rtn;
}

template<class Type>
template<class U>
inline void Cmatrix<Type>::norm_p(Cmatrix<U>&rst,long double p,size_t beg,size_t end,Method flag=LINE) const
{
	if(IsEmpty()) throw bad_oper("norm_p出错:矩阵为空!");
	if(beg==0||end==0) throw bad_oper("norm_p出错:下标不能为0!");
	if(beg>end) std::swap(beg,end);
	if(flag==LINE)
	{
		if(end>line_count) throw bad_oper("norm_p出错:列标出界!");
		rst.resize(1,end-beg+1);
		for(size_t j=beg;j<=end;j++)
		{
			long double rtn=0;
			for(size_t i=1;i<=row_count;i++)
				rtn+=std::pow(long double(aufunc::abs((*this)(i,j))),p);
			rst(1,j-beg+1)=std::pow(rtn,1/p);
		}
	}
	else
	{
		if(end>row_count) throw bad_oper("norm_p出错:行标出界!");
		rst.resize(end-beg+1,1);
		for(size_t j=beg;j<=end;j++)
		{
			long double rtn=0;
			for(size_t i=1;i<=line_count;i++)
				rtn+=std::pow(long double(aufunc::abs((*this)(j,i))),p);
			rst(j-beg+1,1)=std::pow(rtn,1/p);
		}
	}
}

template<class Type>
inline Cmatrix<Type> Cmatrix<Type>::orthoriz(Method flag=LINE) const
{
	Cmatrix<Type> tmp(*this);
	if(!tmp.orthorized(flag)) throw bad_oper("orthoriz出错:向量组线性相关!");
	return tmp;
}

template<class Type>
inline Cmatrix<Type> Cmatrix<Type>::vec_unit(Method flag=LINE) const
{
	Cmatrix<Type> tmp(*this);
	tmp.vec_united(flag);
	return tmp;
}

template<class Type>
inline Cmatrix<Type> Cmatrix<Type>::vec_normaliz(Method flag=LINE) const
{
	Cmatrix<Type> tmp(*this);
	tmp.vec_normalized(flag);
	return tmp;
}

//变动型数学运算************************************************
template<class Type>
inline void Cmatrix<Type>::Msin()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=std::sin((*this)(i,j));
}

template<class Type>
inline void Cmatrix<Type>::Msinh()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=std::sinh((*this)(i,j));
}

template<class Type>
inline void Cmatrix<Type>::Mcos()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=std::cos((*this)(i,j));
}

template<class Type>
inline void Cmatrix<Type>::Mcosh()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=std::cosh((*this)(i,j));
}

template<class Type>
inline void Cmatrix<Type>::Mexp()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=std::exp((*this)(i,j));
}

template<class Type>
inline void Cmatrix<Type>::Mlog()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=std::log((*this)(i,j));
}

template<class Type>
inline void Cmatrix<Type>::Mlog10()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=std::log10((*this)(i,j));
}

template<class Type>
template<class U>
inline void Cmatrix<Type>::Mpow(const U& power)
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=std::pow((*this)(i,j),power);
}

template<class Type>
inline void Cmatrix<Type>::Mabs()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=aufunc::abs((*this)(i,j));
}

template<class Type>
inline void Cmatrix<Type>::Masin()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=std::asin((*this)(i,j));
}

template<class Type>
inline void Cmatrix<Type>::Macos()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=std::acos((*this)(i,j));
}

template<class Type>
inline void Cmatrix<Type>::Matan()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=std::atan((*this)(i,j));
}

template<class Type>
inline void Cmatrix<Type>::Msign()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
		{
			if((*this)(i,j)>0) (*this)(i,j)=1;
			else if((*this)(i,j)==0) (*this)(i,j)=0;
			else (*this)(i,j)=-1;
		}		
}

template<class Type>
inline void Cmatrix<Type>::Mconj()
{
	for (size_t i=1;i<=row_count;++i)
		for (size_t j=1;j<=line_count;++j)
			(*this)(i,j)=aufunc::conj((*this)(i,j))	;	
}