/*
 * seq_tools.hpp
 *
 *  Created on: Jan 25, 2014
 *      Author: martin
 *
 *
 */

#ifndef SEQ_TOOLS_HPP_
#define SEQ_TOOLS_HPP_

#include <armadillo>
#include "samtools/faidx.h"

//TODO print file name on error
std::string read_fasta(
		std::string const& filename,
		t_position offset,
		std::string const& ref,
		t_position pos,
		t_count length) {

	faidx_t * fai = fai_load(filename.c_str());

	if (fai == 0) {
		throw std::runtime_error("read_fasta : unable to open file");
	}

	if(pos < offset) {
		throw std::runtime_error("read_fasta : unable to fetch region");
	}

	int l;
	char * cstr_ref = faidx_fetch_seq(fai, const_cast<char*>(ref.c_str()),
			pos-offset, pos-offset + length-1, &l);

	if (cstr_ref == 0) {

		//cout << ref << " : " << pos-offset << " - " << pos-offset + length-1 << endl;

		throw std::runtime_error("read_fasta : unable to fetch region");
	}

	//Note the first base is at postion 0
	std::string s(cstr_ref);

	//Clean up
	free(cstr_ref);
	fai_destroy(fai);


	return s;

}


double quality_to_epsilon(char quality) {
	return pow(10, -static_cast<int>(quality) / static_cast<double>(10));
}

t_epsilon_quality quality_to_epsilon(t_quality const& quality) {

	t_epsilon_quality epsilon(quality.n_elem);
	t_epsilon_quality::iterator epsilon_itr = epsilon.begin();
		for(t_quality::const_iterator quality_itr = quality.begin(); quality_itr != quality.end(); ++quality_itr, ++epsilon_itr) {
			*epsilon_itr = quality_to_epsilon(*quality_itr);
		}

	return epsilon;
}


char seq_base_to_int(const char & bit4_code) {

	switch (bit4_code) {
		case 1: //A
			return 3;

		case 2: //C
			return 1;

		case 4: //G
			return 2;

		case 8: //T
			return 4;

		case 15: //N
			return 0;

		default:
			//TODO debug guard
            //cout << "seq_base_to_int" << " : " << static_cast<int>(bit4_code) << endl;
			throw std::runtime_error("Error: Unknown 4bit code in bam file");
		}

}

char str_base_conv(const char & str_char) {

	switch (str_char) {
	case 'N':
		return 0;

	case 'C':
		return 1;

	case 'G':
		return 2;

	case 'A':
		return 3;

	case 'T':
		return 4;

	case 'c':
		return 1;

	case 'g':
		return 2;

	case 'a':
		return 3;

	case 't':
		return 4;

	case '-':
		return 0;

	default:
        //TODO openmp save warning
        // std::ostringstream msg;
       // msg << "Ignoring unknown symbol '" << str_char << "'";
       // report_warning(msg.str());

        return 0;
	}
}

t_quality create_quality_vector(std::string const& quality_str) {

	t_quality q(quality_str.size());
	t_quality::iterator q_itr = q.begin();
	for(std::string::const_iterator quality_itr = quality_str.begin(); quality_itr != quality_str.end(); ++quality_itr, ++q_itr) {
		*q_itr = static_cast<int>(*quality_itr);
	}

	return q;
}

t_epsilon_quality create_epsilon_vector(std::string const& quality_str) {

	t_epsilon_quality epsilon(quality_str.size());
	t_epsilon_quality::iterator epsilon_itr = epsilon.begin();
	for(std::string::const_iterator quality_itr = quality_str.begin(); quality_itr != quality_str.end(); ++quality_itr, ++epsilon_itr) {
		*epsilon_itr = quality_to_epsilon(*quality_itr);
	}

	return epsilon;
}

t_seq_bases create_bases_vector(std::string const& bases_str) {
	t_seq_bases bases(bases_str.size());
	t_seq_bases::iterator seq_itr = bases.begin();
	for(std::string::const_iterator bases_itr = bases_str.begin(); bases_itr != bases_str.end(); ++bases_itr, ++seq_itr) {
		*seq_itr = str_base_conv(*bases_itr);
	}

	return bases;
}

/*
  	subsets coding
#   0 Ø
#	1 R
#	2 A1
# 	3 A2
#	4 A3
#	5 A1 A2
#	6 A3 R
#	7 A1 A3
#	8 A2 R
#	9 A1 R
#	10 A2 A3
#	11 A1 A3 R
#	12 A1 A2 A3
#	13 A2 R A3
#	14 A1 R A2
#	15 A2 A1 R A3

	base coding:

	1 = R
	2 = A1
	3 = A2
	4 = A3

	0 = N
*/

char base_union(char subset, char base) {

	if(base == 0) { //check if base = N
		return subset;
	}
	switch (subset) {
	case 0: // empty set
		return base;

	case 1:
		switch (base) {
		case 1:
			return 1;

		case 2:
			return 9;

		case 3:
			return 8;

		case 4:
			return 6;

		default:
			throw std::runtime_error("base_union - Internal error");
		}
		break;

	case 2:
		switch (base) {
		case 1:
			return 9;

		case 2:
			return 2;

		case 3:
			return 5;

		case 4:
			return 7;

		default:
			throw std::runtime_error("base_union - Internal error");
		}
		break;

	case 3:
		switch (base) {
		case 1:
			return 8;

		case 2:
			return 5;

		case 3:
			return 3;

		case 4:
			return 10;

		default:
			throw std::runtime_error("base_union - Internal error");
		}
		break;

	case 4:

		switch (base) {
		case 1:
			return 6;

		case 2:
			return 7;

		case 3:
			return 10;

		case 4:
			return 4;

		default:
			throw std::runtime_error("base_union - Internal error");
		}
		break;

	case 5:

		switch (base) {
		case 1:
			return 14;

		case 2:
			return 5;

		case 3:
			return 5;

		case 4:
			return 12;

		default:
			throw std::runtime_error("base_union - Internal error");
		}
		break;

	case 6:

		switch (base) {
		case 1:
			return 6;

		case 2:
			return 11;

		case 3:
			return 13;

		case 4:
			return 6;

		default:
			throw std::runtime_error("base_union - Internal error");
		}
		break;

	case 7: //G T

		switch (base) {
		case 1: //C
			return 11;

		case 2: //G
			return 7;

		case 3: //A
			return 12;

		case 4: //T
			return 7;

		default:
			throw std::runtime_error("base_union - Internal error");
		}
		break;

	case 8: //A C

		switch (base) {
		case 1: //C
			return 8;

		case 2: //G
			return 14;

		case 3: //A
			return 8;

		case 4: //T
			return 13;

		default:
			throw std::runtime_error("base_union - Internal error");
		}
		break;

	case 9: //G C

		switch (base) {
		case 1: //C
			return 9;

		case 2: //G
			return 9;

		case 3: //A
			return 14;

		case 4: //T
			return 11;

		default:
			throw std::runtime_error("base_union - Internal error");
		}
		break;

	case 10: //A T

		switch (base) {
		case 1: //C
			return 13;

		case 2: //G
			return 12;

		case 3: //A
			return 10;

		case 4: //T
			return 10;

		default:
			throw std::runtime_error("base_union - Internal error");
		}
		break;

	case 11:

		if (base != 3) {
			return 11;
		}

		return 15;

	case 12:

		if (base != 1) {
			return 12;
		}

		return 15;

	case 13:

		if (base != 2) {
			return 13;
		}

		return 15;
	case 14:

		if (base != 4) {
			return 14;
		}

		return 15;

	case 15:
		return 15;

	default:
		throw std::runtime_error("base_union - Internal error");
		break;
	}

}

#endif /* SEQ_TOOLS_HPP_ */
