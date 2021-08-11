#ifndef __Vec__
#define __Vec__

class Vec {
	public:
		Vec();
		Vec(int);
		Vec(int, double);
		Vec(int, double*);
		Vec(const Vec&);
		virtual ~Vec();
		
		//operadores
		const Vec& operator=(const Vec&);
  		Vec operator+(const Vec&)const;
  		Vec operator+=(const Vec&);
  		Vec operator-(const Vec&)const;
  		Vec operator-=(const Vec&);
		double& operator[](int);
		Vec operator-()const;
		Vec operator+()const;
  		Vec operator*(const Vec&) const;
  		Vec operator*(double) const;  
  		

		//metodos
		void Print();
                void SetEntries (int, double*);
		double At(int)const;
		int size()const;
		double dot(const Vec&)const;
		void swap(int, int);

	private:
		int N;
		double *entries;
		
};

#endif
