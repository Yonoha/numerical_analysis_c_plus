// 1-dim non-steady diffusion equation by fem

#include <algorithm> // for std::min
#include <array>     // for std::array
#include <fstream>   // for std::ofstream
#include <iostream>  // for std::endl
#include <sstream>   // for std::string, std::stringstream
#include <vector>    // fos std::vector
#include <boost/assert.hpp> // for BOOST_ASSERT
#include <boost/format.hpp> // for boost::format

namespace fem{
	using myvec = std::vector<double>;

    class FEM final{
		private:
        	enum class boundary_condi_type{
            	DIRICLET = 0,
            	LEFT_NEUMANN = 1,
            	RIGHT_NEUMANN = 2
        	};

            // DIRICLET or LEFT NEUMANN or RIGHT NEUMANN
        	static auto constexpr BCT = boundary_condi_type::LEFT_NEUMANN;

			// boundary condition
        	static auto constexpr D0 = 0.0; // left Dirichlet
    		static auto constexpr D1 = 0.0; // right Direchlet
			static auto constexpr N0 = 1.5; // left Neumann
	    	static auto constexpr N1 = 1.0; // right Neumann

            // coefficient
            static auto constexpr DIFF = 0.01; // diffusion 
            static auto constexpr THETA = 0.5; // theta method

            // for discrete of x and t
 	    	static auto constexpr ELEMENT = 100; // mesh number
        	static auto constexpr LENGTH = 1.0;  // xrange
            static auto constexpr TEND = 10.0;   // trange
    
    	public:
			static auto constexpr TLOOP = 10000;
            static auto constexpr DT = TEND / TLOOP;
    		static auto constexpr TREP = TLOOP / 1000;
			static auto constexpr NODE = ELEMENT + 1;
			static auto constexpr DX = LENGTH / ELEMENT;
			static auto constexpr KAP = DIFF * DT / (DX * DX);
	
    	private:
        	myvec bound_;
        	myvec diag1_;
			myvec diag2_;
	    	myvec left1_;
			myvec left2_;
	    	myvec right1_;
	    	myvec right2_;
			myvec x_;

    	public:
        	FEM()
            	: bound_(NODE, 0.0), diag1_(NODE, 0.0), diag2_(NODE, 0.0), 
                  left1_(NODE, 0.0), left2_(NODE, 0.0), right1_(NODE, 0.0), right2_(NODE, 0.0), x_(NODE, 0.0)
        	{
				for (auto i = 0; i < NODE; i++){
					// initial x
					x_[i] = std::min(DX * static_cast<double>(i), 1.0 - DX * static_cast<double>(i)); 
				}
			}

        	~FEM() = default;

			// copy constructorでcopy禁止
        	FEM(FEM const &dummy) = delete;

			// operator=()でもcopy禁止
        	FEM & operator=(FEM const &dummy) = delete;

			// boundary condition
        	void boundary();
			void boundary2();

			// TDMA method
	    	void tdma();

			// make stiffness matrix
	    	void mat();

			// output for each TLOOP
			bool result_output(int i);
    };
}

int main(){
    // KAP > 0.5のとき強制終了 
    // Neumannの安定性解析より, KAP <= 0.5でないと正しい結果が得られない
    BOOST_ASSERT(fem::FEM::KAP <= 0.5);

    fem::FEM fem_obj;

	fem_obj.result_output(0);

	if (!fem_obj.result_output(0)){
		std::cerr << "output file not open" << std::endl;
        return -1;
	}

	fem_obj.mat();
	fem_obj.boundary();

	for (auto i = 1; i <= fem::FEM::TLOOP; i++){
		fem_obj.boundary2();
		fem_obj.tdma();

		if (!(i % 1000)){
			fem_obj.result_output(i);

			if (!fem_obj.result_output(i)){
				std::cerr << "output file not open" << std::endl;
        		return -1;
			}
		}
	}

	return 0;
}

namespace fem{
	bool FEM::result_output(int i){
		std::stringstream ss;
		std::string name;
		std::ofstream ofs;

		ss << i;
		name = ss.str();
		name = "data_nonste_diff_" + name + ".txt";
		
		ofs.open(name.c_str());
		if (!ofs){
			return false;
		}

		for (auto j = 0; j < NODE; j++){
			ofs << DX * static_cast<double>(j) << " " << x_[j] << std::endl;
		}

		return true;
	}

    void FEM::boundary(){
		switch (BCT){
            case boundary_condi_type::DIRICLET:
            	diag1_[0] = 1.0;
				diag1_[NODE - 1] = 1.0;
                diag2_[0] = 1.0;
                diag2_[NODE - 1] = 1.0;
			    left1_[NODE - 1] = 0.0;
                left2_[NODE - 1] = 0.0;
        	    right1_[0] = 0.0;
                right2_[0] = 0.0;
                break;
		
            case boundary_condi_type::LEFT_NEUMANN:
        	    diag1_[NODE - 1] = 1.0;
                diag2_[NODE - 1] = 1.0;
				left1_[NODE - 1] = 0.0;
                left2_[NODE - 1] = 0.0;
                break;
			
            case boundary_condi_type::RIGHT_NEUMANN:
			    diag1_[0] = 1.0;
                diag2_[0] = 1.0;
        	    right1_[0] = 0.0;
                right2_[0] = 0.0;
                break;
            
            default:
                BOOST_ASSERT(!"error! wrong boundary condition");
                break;
        }
	}

    void FEM::boundary2(){
		switch (BCT){
            case boundary_condi_type::DIRICLET:
                bound_[0] = D0;
        	    bound_[NODE - 1] = D1;
                break;
		
            case boundary_condi_type::LEFT_NEUMANN:
                bound_[0] = diag2_[0] * x_[0] + right2_[0] * x_[1] - N0 * DIFF;
				bound_[NODE - 1] = D1;
                break;
			
            case boundary_condi_type::RIGHT_NEUMANN:
                bound_[0] = D0;  
				bound_[NODE - 1] = left2_[NODE - 1] * x_[NODE - 2] + left2_[NODE - 1] * x_[NODE - 1] + N1 * DIFF;
                break;
            
            default:
                BOOST_ASSERT(!"error! wrong boundary condition");
                break;
        }

		for (auto j = 1; j < NODE - 1; j++){
            bound_[j] = left2_[j] * x_[j - 1] + diag2_[j] * x_[j] + right2_[j] * x_[j + 1];
        }
	}

	void FEM::mat(){
	    for (auto i = 0; i < ELEMENT; i++){
            // temporal 1
            diag1_[i] += DX * 2.0 / 6.0 / DT;
            diag1_[i + 1] += DX * 2.0 / 6.0 / DT;
            left1_[i + 1] += DX * 1.0 / 6.0 / DT;
            right1_[i] += DX * 1.0 / 6.0 / DT;

			// diffusion 1
		    diag1_[i] += THETA * DIFF / DX;
            diag1_[i + 1] += THETA * DIFF / DX; 
		    left1_[i + 1] -= THETA * DIFF / DX;
            right1_[i] -= THETA * DIFF / DX;

            // temporal 2
            diag2_[i] += DX * 2.0 / 6.0 / DT;
            diag2_[i + 1] += DX * 2.0 / 6.0 / DT;
            left2_[i + 1] += DX * 1.0 / 6.0 / DT;
            right2_[i] += DX * 1.0 / 6.0 / DT;

			// diffusion 2
		    diag2_[i] -= (1.0 - THETA) * DIFF / DX;
			diag2_[i + 1] -= (1.0 - THETA) * DIFF / DX;
		    left2_[i + 1] += (1.0 - THETA) * DIFF / DX; 
            right2_[i] += (1.0 - THETA) * DIFF / DX;
        }
	}

	void FEM::tdma(){
		myvec p(NODE, 0.0), q(NODE, 0.0);

		p[0] = -right1_[0] / diag1_[0];
		q[0] = bound_[0] / diag1_[0];

		for (auto i = 1; i < NODE; i++){
		    p[i] = -right1_[i] / (diag1_[i] + left1_[i] * p[i - 1]);
		    q[i] = (bound_[i] - left1_[i] * q[i - 1]) / (diag1_[i] + left1_[i] * p[i - 1]);
	    }

		x_[NODE - 1] = q[NODE - 1];

		for (auto i = NODE - 2; i >= 0; i--){
			x_[i] = p[i] * x_[i + 1] + q[i];
		}
 	}
}
