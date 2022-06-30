#include <fstream>  // for std::ofstream
#include <iostream> // for std::endl
#include <vector>   // fos std::vector
#include <boost/assert.hpp> // for BOOST_ASSERT

namespace fem{
	using myvec = std::vector<double>;

    class FEM final{
    	private:
        	enum class boundary_condi_type{
            	DIRICLET = 0,
            	LEFT_NEUMANN = 1,
            	RIGHT_NEUMANN = 2
        	};
        
        	static auto constexpr BCT = boundary_condi_type::LEFT_NEUMANN;

			// boundary conditions
        	static auto constexpr D0 = 0.0; // left Dirichlet
    		static auto constexpr D1 = 1.0; // right Direchlet
			static auto constexpr N0 = 1.0; // left Neumann
	    	static auto constexpr N1 = 0.5; // right Neumann

	    	static auto constexpr ELEMENT = 100;
        	static auto constexpr LENGTH = 1.0;
	    	static auto constexpr NODE = ELEMENT + 1;
        	static auto constexpr DX = LENGTH / ELEMENT;
	
    	private:
        	myvec bound_;
        	myvec diag_;
	    	myvec left_;
	    	myvec right_;
			myvec x_;

    	public:
        	FEM()
            	: bound_(NODE, 0.0), diag_(NODE, 0.0), left_(NODE, 0.0), right_(NODE, 0.0), x_(NODE, 1.0)
			{}

        	~FEM() = default;

			// コピーコントラクタでcopy禁止
        	FEM(FEM const &dummy) = delete;

			// operator=()でもcopy禁止
        	FEM & operator=(FEM const &dummy) = delete;

			// result output file
			bool result_output();

			// boundary condition
        	void boundary();

			// make stiffness matrix
	    	void mat();;

			// TDMA method
	    	void tdma();
    };
}

int main(){
    fem::FEM fem_obj;

	fem_obj.mat();
	fem_obj.boundary();
	fem_obj.tdma();

	// fem_obj.result_output();
	if (!fem_obj.result_output()){
        std::cerr << "output file not open" << std::endl;
        return -1;
    }

	return 0;
}

namespace fem{
	bool FEM::result_output(){
		std::ofstream ofs("data_Poisson.txt");
    	if (!ofs) {
        	return false;
    	}

		for (auto i = 0; i < NODE; i++){
			ofs << fem::FEM::DX * static_cast<double>(i) << " " << x_[i] << std::endl;
		}

		return true;
	}

    void FEM::boundary(){
		switch (BCT){
            case boundary_condi_type::DIRICLET:
            	diag_[0] = 1.0;
				diag_[NODE - 1] = 1.0;
			    left_[NODE - 1] = 0.0;
        	    right_[0] = 0.0;
				// boundary conditions
				bound_[0] = D0;
        	    bound_[NODE - 1] = D1;
                break;
		
            case boundary_condi_type::LEFT_NEUMANN:
        	    diag_[NODE - 1] = 1.0;
				left_[NODE - 1] = 0.0;
				// boundary conditions
			    bound_[0] -= N0; // left Neumann
				bound_[NODE - 1] = D1;
                break;
			
            case boundary_condi_type::RIGHT_NEUMANN:
			    diag_[0] = 1.0;
        	    right_[0] = 0.0;
				// boundary conditions
        	    bound_[0] = D0;
			    bound_[NODE - 1] -= N1; // right Neumann
                break;
            
            default:
                BOOST_ASSERT(!"error! wrong boundary condition");
                break;
        }
	}

	void FEM::mat(){
	    for (auto i = 0; i < ELEMENT; i++){
			//boundary
		    bound_[i] += (2.0 * x_[i] + 1.0 * x_[i + 1]) * DX / 6.0;
		    bound_[i + 1] += (1.0 * x_[i] + 2.0 * x_[i + 1]) * DX / 6.0;

			// diffusion (DIFF = 1)
		    diag_[i] += 1.0 / DX;
			diag_[i + 1] += 1.0 / DX; 
            right_[i] -= 1.0 / DX;
		    left_[i + 1] -= 1.0 / DX;
        }
	}

	void FEM::tdma(){
		myvec p(NODE, 0.0);
		myvec q(NODE, 0.0);

		p[0] = -right_[0] / diag_[0];
		q[0] = bound_[0] / diag_[0];

		for (auto i = 1; i < NODE; i++){
		    p[i] = -right_[i] / (diag_[i] + left_[i] * p[i - 1]);
		    q[i] = (bound_[i] - left_[i] * q[i - 1]) / (diag_[i] + left_[i] * p[i - 1]);
	    }

		x_[NODE - 1] = q[NODE - 1];

		for (auto i = NODE - 2; i >= 0; i--){
			x_[i] = p[i] * x_[i + 1] + q[i];
		}
 	}
}
