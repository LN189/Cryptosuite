# include<iostream>
#include <string>
#define MAX 10000 // for string
#include <sstream>
#include <cmath>
#include <stdlib.h>
using namespace std;
class BigInteger {
private:
    string number;
    bool sign;
public:
    BigInteger(); // empty constructor initializes zero
    BigInteger(string s); // "string" constructor
    BigInteger(string s, bool sin); // "string" constructor
    BigInteger(int n); // "int" constructor
    BigInteger(long long n);
    void setNumber(string s);
    const string& getNumber(); // retrieves the number
    void setSign(bool s);
    const bool& getSign();
    BigInteger absolute(); // returns the absolute value
    void operator = (BigInteger b);
    bool operator == (BigInteger b);
    bool operator != (BigInteger b);
    bool operator > (BigInteger b);
    bool operator < (BigInteger b);
    bool operator >= (BigInteger b);
    bool operator <= (BigInteger b);
    BigInteger& operator ++(); // prefix
    BigInteger  operator ++(int); // postfix
    BigInteger& operator --(); // prefix
    BigInteger  operator --(int); // postfix
    BigInteger operator + (BigInteger b);
    BigInteger operator - (BigInteger b);
    BigInteger operator * (BigInteger b);
    BigInteger operator / (BigInteger b);
    BigInteger operator % (BigInteger b);
    BigInteger& operator += (BigInteger b);
    BigInteger& operator -= (BigInteger b);
    BigInteger& operator *= (BigInteger b);
    BigInteger& operator /= (BigInteger b);
    BigInteger& operator %= (BigInteger b);
    BigInteger& operator [] (int n);
    BigInteger operator -(); // unary minus sign
    operator string(); // for conversion from BigInteger to string
private:
    bool equals(BigInteger n1, BigInteger n2);
    bool less(BigInteger n1, BigInteger n2);
    bool greater(BigInteger n1, BigInteger n2);
    string add(string number1, string number2);
    string subtract(string number1, string number2);
    string multiply(string n1, string n2);
    pair<string, long long> divide(string n, long long den);
    string toString(long long n);
    long long toInt(string s);
};


class cs{
	public:
	void disp();
	void modadd();
	void modmult();
	void modexpo();
	BigInteger Modexpo(BigInteger a,BigInteger b,BigInteger c);
	void gcd();
	BigInteger Gcd(BigInteger a,BigInteger b);
	void invmod();
	BigInteger Invmod(BigInteger a,BigInteger b);

//	void dislog();
	void Rsa();
//	void fieope();
	};
class rsa
	{public :
	string o;
	long long pr[10],lso,lblock;
	int *am;
	BigInteger *bm;
	BigInteger *block;
	BigInteger *e;
	BigInteger *dd;
	int *bd;
	int *groupd;
	string st;
	int ek;
	long long dk;
	long stepsize;
	long long p,q,n,phy;
	void prime();
	void amtobm();
	void bmtoblock();
	void blocktoe();
	void etodd();
	void ddtobd();
	void bdtogroup();
	rsa()
	{
		ek=3;
		stepsize=0;
		pr[0]=295075153 ;pr[1]=314606869 ;pr[2]=433024223 ;pr[3]=472882027 ;pr[4]=217645199 ;pr[5]=533000389 ;
		pr[6]=633910111 ;pr[7]=797003413 ;pr[8]=899809363 ;pr[9]=961748941 ;
	}
	};
cs c;

void rsa::blocktoe()
{
	BigInteger ebi(ek),nbi(n);
	long long i;
	for(i=0;i<lblock;i++)
	{
		e[i]=c.Modexpo(block[i],ebi,nbi);
	}
}

void rsa::etodd()
{
	BigInteger dbi(dk),nbi(n);
	long long i;
	for(i=0;i<lblock;i++)
	{
		dd[i]=c.Modexpo(e[i],dbi,nbi);
	}
}

void rsa::bdtogroup()
{
	int i,a,m,j;
	for(i=0;i<lso;i++)
	{
		a=0;m=1;
		for(j=0;j<8;j++,m*=2)
		{
			a+=m*bd[i*8+j];
		}
		groupd[i]=a;
	}
}

void rsa::ddtobd()
{
	BigInteger x(0),w(2),a;
	long long i,j=0;
	for(i=0;i<lblock-1;i++)
	{
		a=block[i];
		while(a!=x)
		{
			string strings;
			strings=(a%w).getNumber();
			stringstream convert(strings);
			convert>>bd[i*stepsize+j];
			j++;
			a/=w;
		}
		while(j%stepsize!=0)
		{bd[i*stepsize+j]=0;j++;}
	}
	a=block[i];
	while(a!=x)
	{
		string strings=(a%w).getNumber();
		stringstream convert(strings);
		convert>>bd[i*stepsize+j];
		j++;
		a/=w;
	}
	while(j<lso*8)
	{bd[i*stepsize+j]=0;j++;}
}


void rsa::bmtoblock()
{
	long long k,i,j;
	BigInteger s,x(0),y(1),w(2);
	if((lso*8)%stepsize)
	k=(lblock-1)*stepsize;
	else k=lblock*stepsize;
	for(i=0;i*stepsize<k;i++)
	{
		s=y;block[i]=x;
		for(j=0;j<stepsize;j++)
		{
			block[i]+=(bm[i*stepsize+j]*s);
			s*=w;
		}
	}
	if(k!=lblock*stepsize)
	{
		s=y;block[i]=x;
		for(j=0;i*stepsize+j<8*lso;j++)
		{
			block[i]+=(bm[i*stepsize+j]*s);
			s*=w;
		}
	}
}
			
void rsa::prime()
{
	p=pr[rand()%10];
	q=pr[rand()%10];
}
void rsa::amtobm()
{
	BigInteger x(0),y(1);
	int a,r,k;
	long long i;
	for(i=0;i<8*lso;i+=8)
	{
		k=0;
		a=am[i/8];
		while(a)
		{
			r=a%2;
			if(r)
			bm[i+k]=y;
			else bm[i+k]=x;
			k++;
		}
		while(k<8)
		bm[k++]=x;
	}
}
int main()
{
	cs pro;
	pro.disp();
	return 0;
}


//------------------------------------------------------------------------------

void cs::disp()
{
	string z;
	int choice;
	do
	{
	cout << "\033[2J\033[1;1H";
	cout<<"\t\t\t CS PROJECT \n\n";
	cout<<"\t1.\tModulo Addition\n\t2.\tModulo Multiplication\n\t3.\tModulo exponentiation\n";
	cout<<"\t4.\tGCD\n\t5.\tInverse Modulo\n\t6.\tRSA\n";
	cout<<"\t7.\tExit\n";
	cout<<"\tEnter your choice:";
	cin>>choice; 
	cout << "\033[2J\033[1;1H";
	switch(choice)
	{
		case 1:modadd();break;
		case 2:modmult();break;
		case 3:modexpo();break;
		case 4:gcd();break;
		case 5:invmod();break;
//		case 6:dislog();break;
		case 6:Rsa();break;
//		case 8:fieope();break;
		case 7:break;
		default:cout<<"\t\tInvalid input.Try again";getline(cin.ignore(),z);
	}
	}while(choice!=7);
}
void cs::modmult()
{
	string a,b,c,z;
	cout<<"\t\t\t CS PROJECT \n\n\n";
	cout<<"Enter the first number:";
	cin>>a;
	cout<<"Enter the second number:";
	cin>>b;
	cout<<"Enter the divisor:";
	cin>>c;
	BigInteger ab(a),bb(b),cb(c),d;
	d=(ab*bb)%cb;
	cout<<"The remainder is ="<<d.getNumber();
	getline(cin.ignore(),z);
}

void cs::modadd()
{
	string a,b,c,z;
	cout<<"\t\t\t CS PROJECT \n\n\n";
	cout<<"Enter the first number:";
	cin>>a;
	cout<<"Enter the second number:";
	cin>>b;
	cout<<"Enter the divisor:";
	cin>>c;
	BigInteger ab(a),bb(b),cb(c),d;
	d=((ab%cb)+(bb%cb))%cb;
	cout<<"The remainder is ="<<d.getNumber();
	getline(cin.ignore(),z);
}
void cs::modexpo()
{
	string as,bs,cs,z;
	cout<<"\t\t\t CS PROJECT \n\n\n";
	cout<<"Enter the first number:";
	cin>>as;
	cout<<"Enter the exponential:";
	cin>>bs;
	cout<<"Enter the divisor:";
	cin>>cs;
	BigInteger a(as),b(bs),c(cs),d,x(0);
	d=Modexpo(a,b,c);	
	cout<<"The remainder is ="<<d.getNumber();
	getline(cin.ignore(),z);
}
BigInteger cs::Modexpo(BigInteger a,BigInteger b,BigInteger c)
{
	BigInteger d(1),x(0);
	
	BigInteger radix(2),base(a);
	while(b>x)
	{
		if((b%radix)!=x)
		d*=(b%radix)*(base%c);
		b/=radix;
		base=(base*base)%c;
		d%=c;	
	}
	return d;
}
	
void cs::gcd()
{
	string as,bs,z;
	cout<<"\t\t\tCS PROJECT \n\n\n";
	cout<<"Enter the first number:";
	cin>>as;
	cout<<"Enter the second number:";
	cin>>bs;
	BigInteger a(as),b(bs),c;
	c=Gcd(a,b);
	cout<<"The Gcd is="<<c.getNumber();
	getline(cin.ignore(),z);
}
BigInteger cs::Gcd(BigInteger a,BigInteger b)
{
	BigInteger temp,x(0);
	if(a>b)
	{
		while((a%b)!=x)
		{
			temp=b;
			b=a%b;
			a=temp;
		}
		return b;
	}
	else
	{
		while((b%a)!=x)
		{
			temp=a;
			a=b%a;
			b=temp;
		}
		return a;
	}
}


void cs::invmod()
{
	string as,bs,z;
	cout<<"\t\t\t CS PROJECT \n\n\n";
	cout<<"Enter the first number:";
	cin>>as;
	cout<<"Enter the second number:";
	cin>>bs;
	BigInteger a(as),b(bs),d,x(0);
	d=Invmod(a,b);
	if(d!=x)	
	cout<<"The inverse is ="<<d.getNumber();
	else
	cout<<"Inverse doesn't exist";	
	getline(cin.ignore(),z);
}

BigInteger cs::Invmod(BigInteger a,BigInteger b)
{	BigInteger x(0),y(1);
	if(Gcd(a,b)==y)
	{	
		BigInteger b0 = b, t, q;
		BigInteger x0(0), x1(1),y(1),x(0);
		if (b == y) return y;
		while (a > y) {
		q = a / b;	
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - (q * x0); x1 = t;
		}
		if (x1 < x) x1 += b0;
		return x1;
	}
	else
	return x;
}
void cs::Rsa()
{
	string rs;
	long long iq;
	rsa a;
	a.prime();
	a.n=a.p*a.q;
	a.phy=(a.p-1)*(a.q-1);
	cout<<"\t\t\tCS PROJECT \n\n\n";
	cout<<"Enter the message:";
	cin>>a.o;
	a.lso=a.o.length();
	a.bm=new BigInteger[8*a.lso];
	a.am=new int[8*a.lso];
	long long n=a.n;	
	while(n>0)
	{
		a.stepsize++;
		n/=2;
	}
	if((8*a.lso)%a.stepsize==0)
	a.lblock=(8*a.lso)/a.stepsize;
	else
	a.lblock=((8*a.lso)/a.stepsize)+1;
	a.block=new BigInteger[a.lblock];
	a.e=new BigInteger[a.lblock];
	a.dd=new BigInteger[a.lblock];
	a.bd=new int[8*a.lso];
	a.groupd=new int[a.lso];
	BigInteger eki(a.ek),phyi(a.phy);
	string str=Invmod(eki,phyi).getNumber();
	stringstream convert(str);
	convert>>a.dk;
	a.amtobm();
	a.bmtoblock();
	a.blocktoe();
	cout<<"The encrypted message is:";
	for(iq=0;iq<a.lblock;iq++)
	cout<<a.e[iq].getNumber();
	cout<<endl;
	a.etodd();
	a.ddtobd();
	a.bdtogroup();
	for(iq=0;iq<a.lso;iq++)
	a.st=a.st+char(a.groupd[iq]);
	cout<<"The decrypted message is:"<<a.st;
	getline(cin.ignore(),rs);
}



BigInteger::BigInteger() { // empty constructor initializes zero
    number = "0";
    sign = false;
}

BigInteger::BigInteger(string s) { // "string" constructor
    if( isdigit(s[0]) ) { // if not signed
        setNumber(s);
        sign = false; // +ve
    } else {
        setNumber( s.substr(1) );
        sign = (s[0] == '-');
    }
}

BigInteger::BigInteger(string s, bool sin) { // "string" constructor
    setNumber( s );
    setSign( sin );
}

BigInteger::BigInteger(int n) { // "int" constructor
    stringstream ss;
    string s;
    ss << n;
    ss >> s;


    if( isdigit(s[0]) ) { // if not signed
        setNumber( s );
        setSign( false ); // +ve
    } else {
        setNumber( s.substr(1) );
        setSign( s[0] == '-' );
    }
}
BigInteger::BigInteger(long long n) { // "int" constructor
    stringstream ss;
    string s;
    ss << n;
    ss >> s;


    if( isdigit(s[0]) ) { // if not signed
        setNumber( s );
        setSign( false ); // +ve
    } else {
        setNumber( s.substr(1) );
        setSign( s[0] == '-' );
    }
}

void BigInteger::setNumber(string s) {
    number = s;
}

const string& BigInteger::getNumber() { // retrieves the number
    return number;
}

void BigInteger::setSign(bool s) {
    sign = s;
}

const bool& BigInteger::getSign() {
    return sign;
}

BigInteger BigInteger::absolute() {
    return BigInteger( getNumber() ); // +ve by default
}

void BigInteger::operator = (BigInteger b) {
    setNumber( b.getNumber() );
    setSign( b.getSign() );
}

bool BigInteger::operator == (BigInteger b) {
    return equals((*this) , b);
}

bool BigInteger::operator != (BigInteger b) {
    return ! equals((*this) , b);
}

bool BigInteger::operator > (BigInteger b) {
    return greater((*this) , b);
}

bool BigInteger::operator < (BigInteger b) {
    return less((*this) , b);
}

bool BigInteger::operator >= (BigInteger b) {
    return equals((*this) , b)
           || greater((*this), b);
}

bool BigInteger::operator <= (BigInteger b) {
    return equals((*this) , b)
           || less((*this) , b);
}

BigInteger& BigInteger::operator ++() { // prefix
    (*this) = (*this) + 1;
    return (*this);
}

BigInteger BigInteger::operator ++(int) { // postfix
    BigInteger before = (*this);

    (*this) = (*this) + 1;

    return before;
}

BigInteger& BigInteger::operator --() { // prefix
    (*this) = (*this) - 1;
    return (*this);

}

BigInteger BigInteger::operator --(int) { // postfix
    BigInteger before = (*this);

    (*this) = (*this) - 1;

    return before;
}

BigInteger BigInteger::operator + (BigInteger b) {
    BigInteger addition;
    if( getSign() == b.getSign() ) { // both +ve or -ve
        addition.setNumber( add(getNumber(), b.getNumber() ) );
        addition.setSign( getSign() );
    } else { // sign different
        if( absolute() > b.absolute() ) {
            addition.setNumber( subtract(getNumber(), b.getNumber() ) );
            addition.setSign( getSign() );
        } else {
            addition.setNumber( subtract(b.getNumber(), getNumber() ) );
            addition.setSign( b.getSign() );
        }
    }
    if(addition.getNumber() == "0") // avoid (-0) problem
        addition.setSign(false);

    return addition;
}

BigInteger BigInteger::operator - (BigInteger b) {
    b.setSign( ! b.getSign() ); // x - y = x + (-y)
    return (*this) + b;
}

BigInteger BigInteger::operator * (BigInteger b) {
    BigInteger mul;

    mul.setNumber( multiply(getNumber(), b.getNumber() ) );
    mul.setSign( getSign() != b.getSign() );

    if(mul.getNumber() == "0") // avoid (-0) problem
        mul.setSign(false);

    return mul;
}

// Warning: Denomerator must be within "long long" size not "BigInteger"
BigInteger BigInteger::operator / (BigInteger b) {
    long long den = toInt( b.getNumber() );
    BigInteger div;

    div.setNumber( divide(getNumber(), den).first );
    div.setSign( getSign() != b.getSign() );

    if(div.getNumber() == "0") // avoid (-0) problem
        div.setSign(false);

    return div;
}

// Warning: Denomerator must be within "long long" size not "BigInteger"
BigInteger BigInteger::operator % (BigInteger b) {
    long long den = toInt( b.getNumber() );

    BigInteger rem;
    long long rem_int = divide(number, den).second;
    rem.setNumber( toString(rem_int) );
    rem.setSign( getSign() != b.getSign() );

    if(rem.getNumber() == "0") // avoid (-0) problem
        rem.setSign(false);

    return rem;
}

BigInteger& BigInteger::operator += (BigInteger b) {
    (*this) = (*this) + b;
    return (*this);
}

BigInteger& BigInteger::operator -= (BigInteger b) {
    (*this) = (*this) - b;
    return (*this);
}

BigInteger& BigInteger::operator *= (BigInteger b) {
    (*this) = (*this) * b;
    return (*this);
}

BigInteger& BigInteger::operator /= (BigInteger b) {
    (*this) = (*this) / b;
    return (*this);
}

BigInteger& BigInteger::operator %= (BigInteger b) {
    (*this) = (*this) % b;
    return (*this);
}

BigInteger& BigInteger::operator [] (int n) {
    return *(this + (n*sizeof(BigInteger)));
}

BigInteger BigInteger::operator -() { // unary minus sign
    return (*this) * -1;
}

BigInteger::operator string() { // for conversion from BigInteger to string
    string signedString = ( getSign() ) ? "-" : ""; // if +ve, don't print + sign
    signedString += number;
    return signedString;
}

bool BigInteger::equals(BigInteger n1, BigInteger n2) {
    return n1.getNumber() == n2.getNumber()
           && n1.getSign() == n2.getSign();
}

bool BigInteger::less(BigInteger n1, BigInteger n2) {
    bool sign1 = n1.getSign();
    bool sign2 = n2.getSign();

    if(sign1 && ! sign2) // if n1 is -ve and n2 is +ve
        return true;

    else if(! sign1 && sign2)
        return false;

    else if(! sign1) { // both +ve
        if(n1.getNumber().length() < n2.getNumber().length() )
            return true;
        if(n1.getNumber().length() > n2.getNumber().length() )
            return false;
        return n1.getNumber() < n2.getNumber();
    } else { // both -ve
        if(n1.getNumber().length() > n2.getNumber().length())
            return true;
        if(n1.getNumber().length() < n2.getNumber().length())
            return false;
        return n1.getNumber().compare( n2.getNumber() ) > 0; // greater with -ve sign is LESS
    }
}

bool BigInteger::greater(BigInteger n1, BigInteger n2) {
    return ! equals(n1, n2) && ! less(n1, n2);
}

string BigInteger::add(string number1, string number2) {
    string add = (number1.length() > number2.length()) ?  number1 : number2;
    char carry = '0';
    int differenceInLength = abs( (int) (number1.size() - number2.size()) );

    if(number1.size() > number2.size())
        number2.insert(0, differenceInLength, '0'); // put zeros from left

    else// if(number1.size() < number2.size())
        number1.insert(0, differenceInLength, '0');

    for(int i=number1.size()-1; i>=0; --i) {
        add[i] = ((carry-'0')+(number1[i]-'0')+(number2[i]-'0')) + '0';

        if(i != 0) {
            if(add[i] > '9') {
                add[i] -= 10;
                carry = '1';
            } else
                carry = '0';
        }
    }
    if(add[0] > '9') {
        add[0]-= 10;
        add.insert(0,1,'1');
    }
    return add;
}

string BigInteger::subtract(string number1, string number2) {
    string sub = (number1.length()>number2.length())? number1 : number2;
    int differenceInLength = abs( (int)(number1.size() - number2.size()) );

    if(number1.size() > number2.size())
        number2.insert(0, differenceInLength, '0');

    else
        number1.insert(0, differenceInLength, '0');

    for(int i=number1.length()-1; i>=0; --i) {
        if(number1[i] < number2[i]) {
            number1[i] += 10;
            number1[i-1]--;
        }
        sub[i] = ((number1[i]-'0')-(number2[i]-'0')) + '0';
    }

    while(sub[0]=='0' && sub.length()!=1) // erase leading zeros
        sub.erase(0,1);

    return sub;
}

string BigInteger::multiply(string n1, string n2) {
    if(n1.length() > n2.length())
        n1.swap(n2);

    string res = "0";
    for(int i=n1.length()-1; i>=0; --i) {
        string temp = n2;
        int currentDigit = n1[i]-'0';
        int carry = 0;

        for(int j=temp.length()-1; j>=0; --j) {
            temp[j] = ((temp[j]-'0') * currentDigit) + carry;

            if(temp[j] > 9) {
                carry = (temp[j]/10);
                temp[j] -= (carry*10);
            } else
                carry = 0;

            temp[j] += '0'; // back to string mood
        }

        if(carry > 0)
            temp.insert(0, 1, (carry+'0'));

        temp.append((n1.length()-i-1), '0'); // as like mult by 10, 100, 1000, 10000 and so on

        res = add(res, temp); // O(n)
    }

    while(res[0] == '0' && res.length()!=1) // erase leading zeros
        res.erase(0,1);

    return res;
}

pair<string, long long> BigInteger::divide(string n, long long den) {
    long long rem = 0;
    string result;
    result.resize(MAX);

    for(int indx=0, len = n.length(); indx<len; ++indx) {
        rem = (rem * 10) + (n[indx] - '0');
        result[indx] = rem / den + '0';
        rem %= den;
    }
    result.resize( n.length() );

    while( result[0] == '0' && result.length() != 1)
        result.erase(0,1);

    if(result.length() == 0)
        result = "0";

    return make_pair(result, rem);
}

string BigInteger::toString(long long n) {
    stringstream ss;
    string temp;

    ss << n;
    ss >> temp;

    return temp;
}

long long BigInteger::toInt(string s) {
    long long sum = 0;

    for(int i=0; i<s.length(); i++)
        sum = (sum*10) + (s[i] - '0');

    return sum;
}

