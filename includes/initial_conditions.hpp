#ifndef INITIAL_CONDITIONS
#define INITIAL_CONDITIONS

class initial_conditions{
	private:
		double offset1 = 0.25;
		double offset2 = 0.75;
		double sigma = 0.05;
	public:
		double get_offset1();
		double get_offset2();
		double get_sigma();
		double phi_init(double i_offset1, double i_offset2, double i_sigma, double i_x);
		double psi_init(double i_offset1, double i_offset2, double i_sigma, double i_x);
};

#endif