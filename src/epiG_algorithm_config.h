/* Routines for TODO
 Intended for use with R.
 Copyright (C) 2012 Martin Vincent

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

#ifndef EPIG_ALGORITHM_CONFIG_H_
#define EPIG_ALGORITHM_CONFIG_H_

template<typename type>
static type getConfigAttribute(rList const& config, std::string const& name) {

	int index;
	if (index = config.getIndex(name), index >= 0) {

		return get_value<type>(config.get(index));

	} else {

		std::string msg = "Missing configuration parameter : ";
		throw std::domain_error(msg.append(name).c_str());
		return type(); //avoid compiler warnings
	}
}

//TODO do we need this
template<typename type>
static field<type> getConfigList(rList const& config, std::string const& name) {

	int index;
	if (index = config.getIndex(name), index >= 0) {

		return get_field<type>(config.get(index));

	} else {

		std::string msg = "Missing configuration parameter : ";
		throw std::domain_error(msg.append(name).c_str());
		return field<type>(); //avoid compiler warnings
	}
}

class AlgorithmConfiguration {

public:

	const t_count max_iterations; //max number of iterations per optim loop

	const t_models fwd_model;
	const t_models rev_model;

	const t_models fwd_GpC_model;
	const t_models rev_GpC_model;

    const t_count reads_hard_limit;

	const arma::Col<double> haplochain_log_prior;

    const double ref_prior;

    const t_position min_overlap_length;

    bool const use_paired_reads;

    bool const NOMEseq_mode;

	bool const verbose;

	AlgorithmConfiguration(rList const& config) :

			max_iterations(getConfigAttribute<t_count>(config, "max_iterations")),

			fwd_model(getConfigList<t_model>(config, "fwd_model")),

			rev_model(getConfigList<t_model>(config, "rev_model")),

			fwd_GpC_model(getConfigList<t_model>(config, "fwd_GpC_model")),

			rev_GpC_model(getConfigList<t_model>(config, "rev_GpC_model")),

            reads_hard_limit(getConfigAttribute<t_count>(config, "reads_hard_limit")),

            haplochain_log_prior(getConfigAttribute<arma::Col<double> >(config, "log_haplo_prior")),

            ref_prior(getConfigAttribute<double>(config, "ref_prior")),

            min_overlap_length(getConfigAttribute<t_position>(config, "min_overlap_length")),

			use_paired_reads(getConfigAttribute<bool>(config, "use_paired_reads")),

			NOMEseq_mode(getConfigAttribute<bool>(config, "NOMEseq_mode")),

			verbose(getConfigAttribute<bool>(config, "verbose")) {
	}
};

#endif /* EPIG_ALGORITHM_CONFIG_H_ */
