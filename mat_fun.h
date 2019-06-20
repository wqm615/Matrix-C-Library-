//该文件为Cmatrix类的全局函数

inline long double sin(int a) 
{
	return sin(long double(a));
}

//全局数学函数
template<class Type,class U>
inline Cmatrix<Type> Pmutiply(const Cmatrix<Type>&L,const Cmatrix<U>&R)
{
	if(L.IsEmpty()||R.IsEmpty()) throw bad_oper("Pmutiply出错:矩阵为空!");
	if (L.row()!=R.row()||L.line()!=R.line())
		throw bad_oper("矩阵的行列数不一样！");
	Cmatrix<Type> tmp(L);
	for (size_t i=1;i<=L.row();++i)
		for (size_t j=1;j<=L.line();++j)
			tmp(i,j)*=R(i,j);
	return tmp;
}

template<class Type,class U>
inline Cmatrix<Type> Cmutiply(const Cmatrix<Type>&L,const Cmatrix<U>&R)
{
	if(L.IsEmpty()||R.IsEmpty()) throw bad_oper("Cmutiply出错:矩阵为空!");
	Cmatrix<Type> temp(L.row()*R.row(),L.line()*R.line(),0);
	for (size_t row_l=1;row_l<=L.row();row_l++)
		for (size_t line_l=1;line_l<=L.line();line_l++)
		{
			for (size_t row_r=1;row_r<=R.row();row_r++)
				for (size_t line_r=1;line_r<=R.line();line_r++)
					temp((row_l-1)*R.row()+row_r,(line_l-1)*R.line()+line_r)=L(row_l,line_l)*R(row_r,line_r);
		}
	return temp;
};


template<class Type>
inline Cmatrix<Type > sin(const Cmatrix<Type >&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=sin(temp(row,line));
	return temp;  //复数和其它类型通用

};

template<class Type>
inline Cmatrix<Type > asin(const Cmatrix<Type >&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=asin(temp(row,line));
	return temp;  //复数和其它类型通用

};


template<class Type>
inline Cmatrix<Type> abs(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=aufunc::abs(temp(row,line));
	return temp;
};

template<class Type>
inline Cmatrix<Type> cos(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=cos(temp(row,line));
	return temp;
};

template<class Type>
inline Cmatrix<Type> acos(Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=acos(temp(row,line));
	return temp;
};

template<class Type>
inline Cmatrix<Type> cosh(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=cosh(temp(row,line));
	return temp;
}

template<class Type>
inline Cmatrix<Type> tan(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=tan(temp(row,line));
	return temp;
}

template<class Type>
inline Cmatrix<Type> atan(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=atan(temp(row,line));
	return temp;
}

template<class Type>
inline Cmatrix<Type> tanh(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=tanh(temp(row,line));
	return temp;
}

template<class Type>
inline Cmatrix<Type> sinh(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=sinh(temp(row,line));
	return temp;
}

template<class Type>
inline Cmatrix<Type> exp(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=exp(temp(row,line));
	return temp;
}

template<class Type>
inline Cmatrix<Type> log(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=log(temp(row,line));
	return temp;
}

template<class Type>
inline Cmatrix<Type> log10(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=log10(temp(row,line));
	return temp;
}

template<class Type,class U>
inline Cmatrix<Type> pow(const Cmatrix<Type>&orig,const U &power)
{
	Cmatrix<Type> temp(orig);
	for (int row=1;row<=orig.row();row++)
		for (int line=1;line<=orig.line();line++)
			temp(row,line)=pow(temp(row,line),power);
	return temp;
}

template<class Type>
inline Cmatrix<Type> sign(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	temp.Msign();
	return temp;
}

template<class Type>
inline Cmatrix<Type> conj(const Cmatrix<Type>&orig)
{
	Cmatrix<Type> temp(orig);
	temp.Mconj();
	return temp;
}


//工具函数

 template<class Type,class Type2>
 inline  Cmatrix<Type> operator *(Type2 n, const  Cmatrix<Type> &orig)
		{return orig*n;}

 template<class Type,class Type2>
 inline  Cmatrix<Type> operator +(Type2 n, const  Cmatrix<Type> &orig)
		{return orig+n;}

 template<class Type,class Type2>
 inline  Cmatrix<Type> operator -(Type2 n, const  Cmatrix<Type> &orig)
		{
			Cmatrix<Type> tmp(orig);
			tmp*=-1;
			tmp+=n;
			return tmp;
		}

 template<class Type2,class Type>
 inline  Cmatrix<Type> operator /(Type2 n, const  Cmatrix<Type> &orig)
		{
			Cmatrix<Type> tmp(orig);
			tmp.invered()
			return tmp*=n;
        }

