//const_iterator的实现文件

//构造函数
template<class Type>
inline Cmatrix<Type>::const_iterator::const_iterator()
{ptr=NULL;row=line=0;};

 template<class Type>
 inline Cmatrix<Type>::const_iterator::const_iterator(const typename Cmatrix<Type>::iterator &orig)
		{*this=orig;}

 template<class Type>
 inline Cmatrix<Type>::const_iterator::const_iterator(const const_iterator &orig)
		{*this=orig;}



//赋值运算	
 template<class Type>
 inline typename Cmatrix<Type>::const_iterator & Cmatrix<Type>::const_iterator::operator=(const typename Cmatrix<Type>::iterator &orig)
		{
			ptr=orig.ptr;
			row=orig.row;
			line=orig.line;
			return *this;
		}


 template<class Type>
 inline typename Cmatrix<Type>::const_iterator & Cmatrix<Type>::const_iterator::operator=(const const_iterator &orig)
		{
			ptr=orig.ptr;
			row=orig.row;
			line=orig.line;
			return *this;
		};


 //解引用
 template<class Type>
 inline const Type & Cmatrix<Type>::const_iterator::operator*() const
		{return (*ptr)(row,line);};

 template<class Type>
 inline const Type* Cmatrix<Type>::const_iterator::operator->() const
		{return &((*ptr)(row,line));};

 template<class Type>
 inline const Type& Cmatrix<Type>::const_iterator::operator [](size_t n) const
		{return *((*this)+n);}

 //自加减
 template<class Type>
 inline typename Cmatrix<Type>::const_iterator & Cmatrix<Type>::const_iterator::operator ++()
		{
			if (++row>ptr->row_count)
			{++line;row=1;}
			return *this;
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_iterator Cmatrix<Type>::const_iterator::operator ++(int) 
		{
			const_iterator temp;
			temp=*this;
			if (++row>ptr->row_count)
			{++line;row=1;}
			return temp;
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_iterator & Cmatrix<Type>::const_iterator::operator --()
		{
			if (--row<1)
			{--line;row=ptr->row_count;}
			return *this;
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_iterator Cmatrix<Type>::const_iterator::operator --(int) 
		{
			const_iterator temp;
			temp=*this;
			if (--row<1)
			{--line;row=ptr->row_count;}
			return temp;
		};

 //逻辑运算
 template<class Type>
 inline	bool Cmatrix<Type>::const_iterator::operator ==(const  const_iterator &it) const
		{return  row==it.row&&line==it.line&&ptr==it.ptr;};

 template<class Type>
 inline	bool Cmatrix<Type>::const_iterator::operator ==(const typename Cmatrix<Type>::iterator &it) const
		{return  row==it.row&&line==it.line&&ptr==it.ptr;};

 template<class Type>
 inline bool Cmatrix<Type>::const_iterator::operator<(const  const_iterator &it) const
		{
			if(ptr!=it.ptr) throw bad_oper("迭代器不属于同一个矩阵！");
			return row<it.row||(row==it.row&&line<it.line);
		};

 template<class Type>
 inline bool Cmatrix<Type>::const_iterator::operator<(const typename Cmatrix<Type>::iterator &it) const
		{
			if(ptr!=it.ptr) throw bad_oper("迭代器不属于同一个矩阵！");
			return row<it.row||(row==it.row&&line<it.line);
		};

 template<class Type>
 inline bool Cmatrix<Type>::const_iterator::operator<=(const  const_iterator &it) const
		{return (*this)<it||(*this)==it;};

 template<class Type>
 inline bool Cmatrix<Type>::const_iterator::operator<=(const  typename Cmatrix<Type>::iterator &it) const
		{return (*this)<it||(*this)==it;};

 template<class Type>
 inline bool Cmatrix<Type>::const_iterator::operator>(const  const_iterator &it) const
		{return !((*this)<=it);};

 template<class Type>
 inline bool Cmatrix<Type>::const_iterator::operator>(const  typename Cmatrix<Type>::iterator &it) const
		{return !((*this)<=it);};

 template<class Type>
 inline bool Cmatrix<Type>::const_iterator::operator>=(const const_iterator &it) const
		{return !((*this)<it);};

 template<class Type>
 inline bool Cmatrix<Type>::const_iterator::operator>=(const typename Cmatrix<Type>::iterator &it) const
		{return !((*this)<it);};

 template<class Type>
 inline bool Cmatrix<Type>::const_iterator::operator!=(const const_iterator &it) const
		{return !((*this)==it);};

 template<class Type>
 inline bool Cmatrix<Type>::const_iterator::operator!=(const typename Cmatrix<Type>::iterator &it) const
		{return !((*this)==it);};


 
 //加减运算
 template<class Type>
 inline ptrdiff_t Cmatrix<Type>::const_iterator::operator -(const const_iterator &it2) const
		{
			if(ptr!=it2.ptr) throw bad_oper("迭代器不属于同一个矩阵！");
			return ((line-1)*(ptr->row_count)+row)-((it2.line-1)*(ptr->row_count)+it2.row);
		};

 template<class Type>
 inline ptrdiff_t Cmatrix<Type>::const_iterator::operator -(const typename Cmatrix<Type>::iterator &it2) const
		{
			if(ptr!=it2.ptr) throw bad_oper("迭代器不属于同一个矩阵！");
			return ((line-1)*(ptr->row_count)+row)-((it2.line-1)*(ptr->row_count)+it2.row);
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_iterator Cmatrix<Type>::const_iterator::operator -(difference_type n) const
		{
			const_iterator temp(*this);
			size_t dist=(temp.line-1)*(temp.ptr->row_count)+temp.row-n;
			temp.line=dist/temp.ptr->row_count+1;
			temp.row=dist%temp.ptr->row_count;
			if (temp.row==0)
			{temp.row=temp.ptr->row_count;--temp.line;}
			return temp;
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_iterator Cmatrix<Type>::const_iterator::operator +(difference_type n) const
		{
			const_iterator temp(*this);
			size_t dist=(temp.line-1)*(temp.ptr->row_count)+temp.row+n;
			temp.line=dist/temp.ptr->row_count+1;
			temp.row=dist%temp.ptr->row_count;
			if (temp.row==0)
			{temp.row=temp.ptr->row_count;--temp.line;}
			return temp;
		};


 template<class Type>
 inline typename Cmatrix<Type>::const_iterator & Cmatrix<Type>::const_iterator::operator +=(difference_type n)
		{
			size_t dist=(line-1)*(ptr->row_count)+row+n;
			line=dist/ptr->row_count+1;
			row=dist%ptr->row_count;
			if (row==0)
			{row=ptr->row_count;--line;};
			return *this;
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_iterator & Cmatrix<Type>::const_iterator::operator -=(difference_type n)
		{
			size_t dist=(line-1)*(ptr->row_count)+row-n;
			line=dist/ptr->row_count+1;
			row=dist%ptr->row_count;
			if (row==0)
			{row=ptr->row_count;--line;};
			return *this;
		};

