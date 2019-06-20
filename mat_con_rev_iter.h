//const_reverse_iterator的实现文件

//构造函数
 template<class Type>
 inline Cmatrix<Type>::const_reverse_iterator::const_reverse_iterator()
 {ptr=NULL;row=line=0;};

 template<class Type>
 inline Cmatrix<Type>::const_reverse_iterator::const_reverse_iterator(const const_reverse_iterator &orig)
		{*this=orig;}

 template<class Type>
 inline Cmatrix<Type>::const_reverse_iterator::const_reverse_iterator(const reverse_iterator &orig)
		{*this=orig;}

 template<class Type>
 inline Cmatrix<Type>::const_reverse_iterator::const_reverse_iterator(const typename Cmatrix<Type>::iterator &orig)
		{
			ptr=orig.ptr;
			row=orig.row;
			line=orig.line;
         }

 template<class Type>
 inline Cmatrix<Type>::const_reverse_iterator::const_reverse_iterator(const const_iterator &orig)
		{
			ptr=orig.ptr;
			row=orig.row;
			line=orig.line;
        }

 //赋值运算符的重载
 template<class Type>
 inline typename Cmatrix<Type>::const_reverse_iterator & Cmatrix<Type>::const_reverse_iterator::operator=(const const_reverse_iterator &orig)
		{
			ptr=orig.ptr;
			row=orig.row;
			line=orig.line;
			return *this;
		}

 template<class Type>
 inline typename Cmatrix<Type>::const_reverse_iterator & Cmatrix<Type>::const_reverse_iterator::operator=(const reverse_iterator &orig)
		{
			ptr=orig.ptr;
			row=orig.row;
			line=orig.line;
			return *this;
		}

 //解引用的重载
 template<class Type>
 inline const Type & Cmatrix<Type>::const_reverse_iterator::operator*() const
		{	
			size_t line_temp,row_temp;
			line_temp=line,row_temp=row-1;
			if (row_temp==0)
			{
				row_temp=ptr->row_count;
				--line_temp;
			}
			return (*ptr)(row_temp,line_temp);
		};

 template<class Type>
 inline const Type* Cmatrix<Type>::const_reverse_iterator::operator->() const
		{
			size_t line_temp,row_temp;
			line_temp=line,row_temp=row-1;
			if (row_temp==0)
			{
				row_temp=ptr->row_count;
				--line_temp;
			}
			return &(*ptr)(row_temp,line_temp);
         };

  template<class Type>
  inline const Type& Cmatrix<Type>::const_reverse_iterator::operator [](size_t n) const
		{
			return *((*this)+n);
		};

 //自加减运算符的重载
 template<class Type>
 inline typename Cmatrix<Type>::const_reverse_iterator&  Cmatrix<Type>::const_reverse_iterator::operator ++()
		{
			if (--row==0)
			{--line;row=ptr->row_count;}
			return *this;
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_reverse_iterator Cmatrix<Type>::const_reverse_iterator::operator ++(int)
		{
			const_reverse_iterator temp;
			temp=*this;
			if (--row==0)
			{--line;row=ptr->row_count;}
			return temp;
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_reverse_iterator & Cmatrix<Type>::const_reverse_iterator::operator --()
		{
			if (++row>ptr->row_count)
			{++line;row=1;}
			return *this;
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_reverse_iterator Cmatrix<Type>::const_reverse_iterator::operator --(int)
		{
			const_reverse_iterator temp;
			temp=*this;
			if (++row>ptr->row_count)
			{++line;row=1;}
			return temp;
		};

 //逻辑运算符的重载
template<class Type>
 inline	bool Cmatrix<Type>::const_reverse_iterator::operator ==(const reverse_iterator &it) const
		{return  row==it.row&&line==it.line&&ptr==it.ptr;};

 template<class Type>
 inline	bool Cmatrix<Type>::const_reverse_iterator::operator ==(const const_reverse_iterator &it) const
		{return  row==it.row&&line==it.line&&ptr==it.ptr;};

 template<class Type>
 inline bool Cmatrix<Type>::const_reverse_iterator::operator<(const reverse_iterator &it) const
		{
			if(ptr!=it.ptr) throw bad_oper("迭代器不属于同一个矩阵！");
			return row>it.row||(row==it.row&&line>it.line);
		};

 template<class Type>
 inline bool Cmatrix<Type>::const_reverse_iterator::operator<(const const_reverse_iterator &it) const
		{
			if(ptr!=it.ptr) throw bad_oper("迭代器不属于同一个矩阵！");
			return row>it.row||(row==it.row&&line>it.line);
		};


 template<class Type>
 inline bool Cmatrix<Type>::const_reverse_iterator::operator<=(const reverse_iterator &it) const
		{return (*this)<it||(*this)==it;};

 template<class Type>
 inline bool Cmatrix<Type>::const_reverse_iterator::operator<=(const const_reverse_iterator &it) const
		{return (*this)<it||(*this)==it;};

 template<class Type>
 inline bool Cmatrix<Type>::const_reverse_iterator::operator>(const reverse_iterator &it) const
		{return !((*this)<=it);};

 template<class Type>
 inline bool Cmatrix<Type>::const_reverse_iterator::operator>(const const_reverse_iterator &it) const
		{return !((*this)<=it);};

 template<class Type>
 inline bool Cmatrix<Type>::const_reverse_iterator::operator>=(const reverse_iterator &it) const
		{return !((*this)<it);};

 template<class Type>
 inline bool Cmatrix<Type>::const_reverse_iterator::operator>=(const const_reverse_iterator &it) const
		{return !((*this)<it);};

 template<class Type>
 inline bool Cmatrix<Type>::const_reverse_iterator::operator!=(const reverse_iterator &it) const
		{return !((*this)==it);};

 template<class Type>
 inline bool Cmatrix<Type>::const_reverse_iterator::operator!=(const const_reverse_iterator &it) const
		{return !((*this)==it);};

 //加减运算符的重载
 template<class Type>
 inline ptrdiff_t  Cmatrix<Type>::const_reverse_iterator::operator -(const const_reverse_iterator &it2) const
		{
			if(ptr!=it2.ptr) throw bad_oper("迭代器不属于同一个矩阵！");
			return ((it2.line-1)*(ptr->row_count)+it2.row)-((line-1)*(ptr->row_count)+row);
		};

 template<class Type>
 inline ptrdiff_t  Cmatrix<Type>::const_reverse_iterator::operator -(const reverse_iterator &it2) const
		{
			if(ptr!=it2.ptr) throw bad_oper("迭代器不属于同一个矩阵！");
			return ((it2.line-1)*(ptr->row_count)+it2.row)-((line-1)*(ptr->row_count)+row);
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_reverse_iterator  Cmatrix<Type>::const_reverse_iterator::operator -(difference_type n) const
		{
			const_reverse_iterator temp(*this);
			size_t dist=(temp.line-1)*(temp.ptr->row_count)+temp.row+n;
			temp.line=dist/temp.ptr->row_count+1;
			temp.row=dist%temp.ptr->row_count;
			if (temp.row==0)
			{temp.row=temp.ptr->row_count;--temp.line;}
			return temp;
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_reverse_iterator  Cmatrix<Type>::const_reverse_iterator::operator +(difference_type n) const
		{
			const_reverse_iterator temp(*this);
			size_t dist=(temp.line-1)*(temp.ptr->row_count)+temp.row-n;
			temp.line=dist/temp.ptr->row_count+1;
			temp.row=dist%temp.ptr->row_count;
			if (temp.row==0)
			{temp.row=temp.ptr->row_count;--temp.line;}
			return temp;
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_reverse_iterator & Cmatrix<Type>::const_reverse_iterator::operator +=(difference_type n)
		{
			size_t dist=(line-1)*(ptr->row_count)+row-n;
			line=dist/ptr->row_count+1;
			row=dist%ptr->row_count;
			if (row==0)
			{row=ptr->row_count;--line;};
			return *this;
		};

 template<class Type>
 inline typename Cmatrix<Type>::const_reverse_iterator & Cmatrix<Type>::const_reverse_iterator::operator -=(difference_type n)
		{
			size_t dist=(line-1)*(ptr->row_count)+row+n;
			line=dist/ptr->row_count+1;
			row=dist%ptr->row_count;
			if (row==0)
			{row=ptr->row_count;--line;};
			return *this;
		};

 //返回正常迭代器
 template<class Type>
 inline typename Cmatrix<Type>::const_iterator Cmatrix<Type>::const_reverse_iterator::base() const
 {
	typename Cmatrix<Type>::const_iterator rtn;
	 rtn.ptr=ptr;
	 rtn.line=line;
	 rtn.row=row;
	 return rtn;
 }