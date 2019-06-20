//该文件为Cmatrix类的定义文件

class bad_oper:public std::logic_error				//异常类
{
public:bad_oper(const std::string & str):std::logic_error(str){};
};

template<class T>									//比较准则
class compare_norm:public std::binary_function<T,T,bool>
{
public:
	bool operator()(const T & elem1,const T &elem2)
	{
		return aufunc::abs(elem1)<aufunc::abs(elem2);
	};
};

enum Method{UNIT,RAND,RANDSYM,RANDHEMIT,LINE,ROW,INC,DEC,EQ,NE,LT,GT,LE,GE,FW,BW,M1,F,G,M,LZY,OTHER};

template<class Type>
class Cmatrix;

//在这添加绑定的友元模板函数的声明
  /*示例:
  template<class Type>
  Cmatrix<Type> Trans(const Cmatrix<Type>&);
  */

template<class Type>
class Cmatrix
{
public://类型定义
	typedef Type value_type;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;
	typedef Type & reference;
	typedef const Type & const_reference;
	typedef Type * pointer;
	typedef const Type * const_pointer;

public://视为0的无限小数
	static long double precision;

private:
	std::vector<std::vector<Type> > matrix;
	size_type row_count;
	size_type line_count;

private:
	inline void _rand(Type v1,Type v2,Method=OTHER);
	inline void _randc(Type v1,Type v2,Method=OTHER);
	inline void _randcu(Type v1,Type v2,Method=OTHER);

public://构造函数
	explicit inline Cmatrix(size_type row=0,size_type line=0,Type elem=Type());
	inline Cmatrix(Method,size_type=1);
	inline Cmatrix(Method,Type,Type,size_t=1,size_t=1);
	template<class U>
	inline Cmatrix(const U* ptr,size_type num,size_type row,size_type line,Type elem=Type(),Method=LINE);
	

	template<class U>
	inline Cmatrix(const Cmatrix<U> &orig);

	template<class U>
	inline Cmatrix(const Cmatrix<U> &orig,size_t pos1,size_t pos2,Method=LINE);

	template<class U>
	inline Cmatrix(const Cmatrix<U> &orig,size_t row1,size_t line1,size_t row2,size_t line2);
	
	template<class InputIterator>
    inline Cmatrix(InputIterator beg,InputIterator end,size_type row,size_type line,const Type &elem=Type(),Method=LINE);


public://非变动型查询
	inline size_type row() const;
	inline size_type line() const;
	inline const Type & operator ()(size_type,size_type) const;
	inline std::pair<size_type,size_type> size() const;
	inline Cmatrix<Type> getrow(size_t beg=1,size_t end=1) const;
	inline Cmatrix<Type> getline(size_t beg=1,size_t end=1) const;
	inline Cmatrix<Type> getzone(size_t row1,size_t line1,size_t row2,size_t line2) const;

	inline void Lmax(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg=1,size_t end=1) const;
	inline void Lmaxa(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg=1,size_t end=1) const;
	inline void Rmax(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg=1,size_t end=1) const;
	inline void Rmaxa(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg=1,size_t end=1) const;
	inline void Mmax(Type& rst,std::pair<size_t,size_t>& pos,Method=LINE) const;
	inline void Mmaxa(Type& rst,std::pair<size_t,size_t>& pos,Method=LINE) const;
	inline void Lmin(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg=1,size_t end=1) const;
	inline void Lmina(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg=1,size_t end=1) const;
	inline void Rmin(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg=1,size_t end=1) const;
	inline void Rmina(Cmatrix<Type>& rst,Cmatrix<size_t>& pos,size_t beg=1,size_t end=1) const;
	inline void Mmin(Type& rst,std::pair<size_t,size_t>& pos,Method=LINE) const;
	inline void Mmina(Type& rst,std::pair<size_t,size_t>& pos,Method=LINE) const;

	inline std::pair<Type,size_t> Lmax(const size_t pos,size_t beg=1,size_t end=1) const;
	inline std::pair<Type,size_t> Lmaxa(const size_t pos,size_t beg=1,size_t end=1) const;
	inline std::pair<Type,size_t> Rmax(const size_t pos,size_t beg=1,size_t end=1) const;
	inline std::pair<Type,size_t> Rmaxa(const size_t pos,size_t beg=1,size_t end=1) const;
	inline std::pair<Type,size_t> Lmin(const size_t pos,size_t beg=1,size_t end=1) const;
	inline std::pair<Type,size_t> Lmina(const size_t pos,size_t beg=1,size_t end=1) const;
	inline std::pair<Type,size_t> Rmin(const size_t _pos,size_t beg=1,size_t end=1) const;
	inline std::pair<Type,size_t> Rmina(const size_t _pos,size_t beg=1,size_t end=1) const;

	inline void Lfind(Cmatrix<size_t>&rst,Type v1,Type v2,std::pair<size_t,size_t>lpos=std::pair<size_t,size_t>(1,1),std::pair<size_t,size_t>rpos=std::pair<size_t,size_t>(1,1),Method=FW) const;
	inline void Lfind(Cmatrix<size_t>&rst,Type v1,Method op,std::pair<size_t,size_t>lpos=std::pair<size_t,size_t>(1,1),std::pair<size_t,size_t>rpos=std::pair<size_t,size_t>(1,1),Method=FW) const;
	inline void Lfinda(Cmatrix<size_t>&rst,long double v1,long double v2,std::pair<size_t,size_t>lpos=std::pair<size_t,size_t>(1,1),std::pair<size_t,size_t>rpos=std::pair<size_t,size_t>(1,1),Method=FW) const;
	inline void Lfinda(Cmatrix<size_t>&rst,double v1,Method op,std::pair<size_t,size_t>lpos=std::pair<size_t,size_t>(1,1),std::pair<size_t,size_t>rpos=std::pair<size_t,size_t>(1,1),Method=FW) const;

	inline void Rfind(Cmatrix<size_t>&rst,Type v1,Type v2,std::pair<size_t,size_t>rpos=std::pair<size_t,size_t>(1,1),std::pair<size_t,size_t>lpos=std::pair<size_t,size_t>(1,1),Method=FW) const;
	inline void Rfind(Cmatrix<size_t>&rst,Type v1,Method op,std::pair<size_t,size_t>rpos=std::pair<size_t,size_t>(1,1),std::pair<size_t,size_t>lpos=std::pair<size_t,size_t>(1,1),Method=FW) const;
	inline void Rfinda(Cmatrix<size_t>&rst,long double v1,long double v2,std::pair<size_t,size_t>rpos=std::pair<size_t,size_t>(1,1),std::pair<size_t,size_t>lpos=std::pair<size_t,size_t>(1,1),Method=FW) const;
	inline void Rfinda(Cmatrix<size_t>&rst,long double v1,Method op,std::pair<size_t,size_t>rpos=std::pair<size_t,size_t>(1,1),std::pair<size_t,size_t>lpos=std::pair<size_t,size_t>(1,1),Method=FW) const;

	inline std::pair<size_t,size_t> Mfind(Type v1,Type v2,Method=LINE,Method=FW) const;
	inline std::pair<size_t,size_t> Mfind(Type value,Method=EQ,Method=LINE,Method=FW) const;
	inline std::pair<size_t,size_t> Mfinda(long double v1,long double v2,Method=LINE,Method=FW) const;
	inline std::pair<size_t,size_t> Mfinda(long double value,Method=EQ,Method=LINE,Method=FW) const;

	template<class InputIterator,class U>
	inline Type interp(InputIterator lit,InputIterator rit,const std::pair<U,U> cord) const;
	template<class InputIterator,class U>
	inline void interp(InputIterator lit,InputIterator rit,const Cmatrix<std::pair<U,U> >&cord,Cmatrix<Type>& rst) const;
	template<class U>
	inline Type interp(const Cmatrix<std::pair<U,U> > &cord,const std::pair<U,U>pos) const;
	template<class U>
	inline void interp(const Cmatrix<std::pair<U,U> > &cord,const Cmatrix<std::pair<U,U> >pos,Cmatrix<Type>&rst) const;

	template<class U>
	inline void real(Cmatrix<U>&rst) const;

	inline Cmatrix<long double> real() const;

	template<class U>
	inline void imag(Cmatrix<U>&rst) const;

	inline Cmatrix<long double> imag() const;

	inline bool IsEmpty() const;

public://变动型矩阵变换
	inline Type & operator ()(size_type,size_type);
	inline void  resize(const size_type row,const size_type line,Type=Type());
	inline void  Lresize(size_type line,Type=Type());
	inline void  Rresize(size_type row,Type=Type());
	template<class U>
	inline void combine(const Cmatrix<U>&orig,Method=LINE);
	inline void row(size_t beg,size_t end);
	inline void line(size_t beg,size_t end);
	inline void zone(size_t row1,size_t line1,size_t row2,size_t line2);

	template<class InputIterator>
	inline void Linsert(const size_t pos,const size_t num,InputIterator beg,InputIterator end,const Type &elem=Type());
	template<class InputIterator>
	inline void Linsert(const size_t pos,const size_t num,InputIterator beg);
	template<class U>
	inline void Linsert(const size_t pos,const size_t line_num,const U* ptr,size_t num,const Type&elem=Type());
    
	template<class InputIterator>
	inline void Rinsert(const size_t pos,const size_t num,InputIterator beg,InputIterator end,const Type&elem=Type());
	template<class InputIterator>
	inline void Rinsert(const size_t pos,const size_t num,InputIterator beg);
	template<class U>
	inline void Rinsert(const size_t pos,const size_t line_num,const U* ptr,size_t num,const Type&elem=Type());

	inline void Lerase(const size_t pos);
	inline void Lerase(size_t beg,size_t end);
	inline void Rerase(const size_t pos);
	inline void Rerase(size_t beg,size_t end);
	inline void Radd(const size_t pos1,const size_t pos2,const size_t pos_r,const Type factor1=Type(1),const Type factor2=Type(1));
	inline void Ladd(const size_t pos1,const size_t pos2,const size_t pos_r,const Type factor1=Type(1),const Type factor2=Type(1));

	template<class InputIterator>
	inline void Lpush_back(InputIterator beg); 
	template<class InputIterator>
	inline void Lpush_back(InputIterator beg,InputIterator end,const Type& elem=Type()); 
	template<class U>
	inline void Lpush_back(const U* ptr,size_t num,const Type& elem=Type()); 

	template<class InputIterator>
	inline void Rpush_back(InputIterator beg); 
	template<class InputIterator>
	inline void Rpush_back(InputIterator beg,InputIterator end,const Type&elem=Type()); 
	template<class U>
	inline void Rpush_back(const U* ptr,size_t num,const Type&elem=Type()); 

	inline void pop_back(Method=LINE);

	template<class InputIterator>
	inline void Lpush_up(InputIterator beg); 
	template<class InputIterator>
	inline void Lpush_up(InputIterator beg,InputIterator end,const Type& elem=Type()); 
	template<class U>
	inline void Lpush_up(const U* ptr,size_t num,const Type& elem=Type()); 

	template<class InputIterator>
	inline void Rpush_up(InputIterator beg); 
	template<class InputIterator>
	inline void Rpush_up(InputIterator beg,InputIterator end,const Type& elem=Type()); 
	template<class U>
	inline void Rpush_up(const U* ptr,size_t num,const Type& elem=Type()); 

	inline void pop_up(Method=LINE);

	inline void Mswap(Cmatrix<Type> &orig);
    inline void Rswap(const size_t pos1,const size_t pos2);
	inline void Lswap(const size_t pos1,const size_t pos2);
	inline void Mscale(const Type factor,bool flag=true);
	inline void Rscale(const size_t pos,const Type factor,bool=true);
	inline void Lscale(const size_t pos,const Type factor,bool=true);
	inline void Rsort(size_t beg=1,size_t end=1,Method=INC);
	inline void Lsort(size_t beg=1,size_t end=1,Method=INC);
	inline void Msort(Method=LINE,Method=INC);
	inline void Rsorta(size_t beg=1,size_t end=1,Method=INC);
	inline void Lsorta(size_t beg=1,size_t end=1,Method=INC);
	inline void Msorta(Method=LINE,Method=INC);
	inline void sort(const size_t *,size_t,Method=INC,Method=LINE);
	inline void sorta(const size_t *,size_t,Method=INC,Method=LINE);

	inline void rand(Type v1,Type v2,Method=OTHER);

	inline void clear();
	inline void Lclear();
	inline void Rclear();


public://赋值操作
	template<class U>
	inline void assign(const U *ptr,size_type num,Method=LINE);
	template <class InputIterator>
	inline void assign(InputIterator beg,InputIterator end,Method=LINE);
	template <class InputIterator>
	inline void assign(InputIterator beg,Method=LINE);
	template<class U>
	inline void assign(size_t row1,size_t line1,size_t row2,size_t line2,const Cmatrix<U>&orig,const size_t row=1,const size_t line=1);
	template<class U>
	inline Cmatrix<Type> & operator =(const Cmatrix<U> &orig);
	inline Cmatrix<Type> & operator =(const Type val);
	
public://数学运算符的重载1
	template<class U>
    inline Cmatrix<Type> & operator +=(const Cmatrix<U> &orig);
	template<class U>
	inline Cmatrix<Type> & operator -=(const Cmatrix<U> &orig);
	template<class U>
	inline Cmatrix<Type> & operator *=(const Cmatrix<U> &orig);
	template<class U>
	inline Cmatrix<Type> & operator /=(const Cmatrix<U> &orig);
	template<class U>
    inline Cmatrix<Type>  operator +(const Cmatrix<U> &orig) const;
	template<class U>
	inline Cmatrix<Type>  operator -(const Cmatrix<U> &orig) const;
	template<class U>
	inline Cmatrix<Type>  operator *(const Cmatrix<U> &orig) const;
	template<class U>
	inline Cmatrix<Type>  operator /(const Cmatrix<U> &orig) const;


public://数学运算符的重载2
    inline  Cmatrix<Type> & operator +=(const Type val);
	inline Cmatrix<Type> & operator -=(const Type val);
	inline Cmatrix<Type> & operator *=(const Type val);
	inline Cmatrix<Type> & operator /=(const Type val);
    inline Cmatrix<Type>  operator +(const Type val) const;
	inline Cmatrix<Type>  operator -(const Type val) const;
	inline Cmatrix<Type>  operator *(const Type val) const;
	inline Cmatrix<Type>  operator /(const Type val) const;


public://输入输出
	template<class elem,class traits>
    inline friend std::basic_istream<elem, traits>&operator>> (std::basic_istream<elem, traits>& is, Cmatrix<Type>&orig);
	template<class elem,class traits>
	inline friend std::basic_ostream<elem, traits>&operator<< (std::basic_ostream<elem, traits>& os, const Cmatrix<Type>&orig);
	
public://变动型矩阵运算

	template<class U>
	inline void add(size_t row1,size_t line1,size_t row2,size_t line2,const Cmatrix<U>&orig,const size_t row=1,const size_t line=1);

	template<class U>
	inline void sub(size_t row1,size_t line1,size_t row2,size_t line2,const Cmatrix<U>&orig,const size_t row=1,const size_t line=1);

	template<class U>
	inline bool Pmutiply(const Cmatrix<U> & orig);

	template<class U>
	inline void LCmutiply(const Cmatrix<U> & orig);

	template<class U>
	inline void RCmutiply(const Cmatrix<U> & orig);

	template<class U>
	inline bool Lmutiply(const Cmatrix<U> & orig);

	template<class U>
	inline bool Rmutiply(const Cmatrix<U> & orig);

	template<class U,class V>
	inline bool mutiply(const Cmatrix<U> &L,const Cmatrix<V> &R);

	inline void united(size_t size);

	inline void transed();

	inline void hermited();

	inline bool invered();

	inline void revered(Method=LINE);

	inline void stded();

	inline void gaussed();

	inline void zeroed();

	inline bool orthorized(Method=LINE);

	inline void vec_united(Method=LINE);

	inline void vec_normalized(Method=LINE);

	inline bool eig_jacob(Cmatrix<Type>& eig,Cmatrix<Type>& vec,Method=INC);  //不能用于复数且必须为对称矩阵
	inline bool eig_subspace(Cmatrix<Type>& eig,Cmatrix<Type>& vec,const size_t num=1); //不能用于复数且必须为对称矩阵

	template<class U>
	inline int solve(const Cmatrix<U>& rval,Cmatrix<Type>& rst);

	template<class InputIterator>
	inline int solve(InputIterator beg,Cmatrix<Type>& rst);
	
	template<class U>
	inline int solve(const U * ptr,size_t num,Cmatrix<Type>& rst);
	
	template<class InputIterator>
	inline int solve(InputIterator beg ,InputIterator end,Cmatrix<Type>& rst);

	inline bool solve_g(const Cmatrix<Type>& rval,Cmatrix<Type>& rst);
	inline bool solve_dolt(const Cmatrix<Type>& rval,Cmatrix<Type>& rst);
	inline bool solve_chlsky(const Cmatrix<Type>& rval,Cmatrix<Type>& rst);
	inline bool solve_e(const Cmatrix<Type>& rval,Cmatrix<Type>& rst); //如果单元对称则使用chlsky求解,否则使用dolt求解


public://非变动型矩阵运算

    inline bool issym() const;
    inline bool ishmt() const;
    inline size_t rank() const;
    inline Cmatrix<Type> trans() const;
    inline Cmatrix<Type> inver() const;
    inline Cmatrix<Type> hermit() const;
    inline Cmatrix<Type> rever(Method=LINE) const;
	inline Cmatrix<Type>  std() const;
	inline void fft(Cmatrix<std::complex<Type> >& rst,Method=LINE) const;
	inline void fft(Cmatrix<Type>& rst,Method=LINE) const;
	inline void fftn(Cmatrix<std::complex<Type> >& rst,Method=LINE) const;
	inline void fftn(Cmatrix<Type>& rst,Method=LINE) const;
	inline Type mode() const;
	inline int quadform() const;
	inline Cmatrix<Type> gauss() const;
	inline bool dolt(Cmatrix<Type>& L,Cmatrix<Type>& R,Method=OTHER) const;
	inline bool crout(Cmatrix<Type>& L,Cmatrix<Type>& R,Method=OTHER) const;
	inline bool chlsky(Cmatrix<Type>& L,Cmatrix<Type>& D) const;
	inline Cmatrix<Type> zero() const;
	inline long double radius(Type=Type(),const size_t=100) const; //不能用于有两个相等最大特征值的矩阵
	inline bool eig_jacob(Cmatrix<Type>& eig,Cmatrix<Type>& vec,Method=INC) const;  //不能用于复数且必须为对称矩阵
	inline bool eig_subspace(Cmatrix<Type>& eig,Cmatrix<Type>& vec,const size_t num=1) const; //不能用于复数且必须为对称矩阵

	template<class U>
	inline int solve(const Cmatrix<U>& rval,Cmatrix<Type>& rst) const;

	template<class InputIterator>
	inline int solve(InputIterator beg,Cmatrix<Type>& rst) const;
	
	template<class U>
	inline int solve(const U * ptr,size_t num,Cmatrix<Type>& rst) const;
	
	template<class InputIterator>
	inline int solve(InputIterator beg ,InputIterator end,Cmatrix<Type>& rst) const;

	inline bool solve_g(const Cmatrix<Type>& rval,Cmatrix<Type>& rst) const;
	inline bool solve_dolt(const Cmatrix<Type>& rval,Cmatrix<Type>& rst) const;
	inline bool solve_chlsky(const Cmatrix<Type>& rval,Cmatrix<Type>& rst) const;
	inline bool solve_e(const Cmatrix<Type>& rval,Cmatrix<Type>& rst) const; //如果单元对称则使用chlsky求解,否则使用dolt求解

	inline long double norm(Method) const;
	inline long double norm(int=1) const;
	template<class U>
	inline void norm_p(Cmatrix<U>&rst,long double p=1,size_t beg=1,size_t end=1,Method=LINE) const;

	inline Cmatrix<Type> orthoriz(Method=LINE) const; //不能用于复数

	inline Cmatrix<Type> vec_unit(Method=LINE) const;

	inline Cmatrix<Type> vec_normaliz(Method=LINE) const;
    
public://变动型数学运算
	inline void Msin();
	inline void Msinh();
	inline void Mcos();
	inline void Mcosh();
	inline void Mexp();
	inline void Mlog();
	inline void Mlog10();
	template<class U>
	inline void Mpow(const U&);
	inline void Mabs();
	inline void Masin();
	inline void Macos();
	inline void Matan();
	inline void Msign();
	inline void Mconj();

};






