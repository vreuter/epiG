#ifndef EPIG_FETCH_R_INTERFACE_HPP_
#define EPIG_FETCH_R_INTERFACE_HPP_

extern "C" {

SEXP r_epiG_fetch_reads(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end);
SEXP r_epiG_fetch_reads_raw(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end);
SEXP r_epiG_fetch_reads_info(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end);

SEXP r_epiG_fetch_header(SEXP r_filename);

SEXP r_epiG_read_fasta(SEXP r_filename, SEXP r_ref, SEXP r_position, SEXP r_length);

}


SEXP epiG_read_fasta(SEXP r_filename, SEXP r_ref, SEXP r_position, SEXP r_length) {

	const std::string filename = get_value<std::string>(r_filename);
	const std::string ref = get_value<std::string>(r_ref);
	const int pos = get_value<int>(r_position);
	const int length = get_value<int>(r_length);

	//Note the first base is at postion 0
	return rObject(create_bases_vector(read_fasta(filename, ref, pos, length)));
}


SEXP r_epiG_read_fasta(SEXP r_filename, SEXP r_ref, SEXP r_position, SEXP r_length) {

	try {

		return epiG_read_fasta(r_filename, r_ref, r_position, r_length);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if (e.what() != NULL) {
			report_error(e.what());
		}

		else {
			report_error("Unknown error");
		}

	} catch (...) {
		report_error("Unknown error");
	}

	return R_NilValue; //Avoid compiler warnings
}



SEXP epiG_fetch_reads(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end) {

	const std::string filename = get_value<std::string>(r_filename);
	const std::string refName = get_value<std::string>(r_refName);
	const t_position start = get_value<t_position>(r_start);
	const t_position end = get_value<t_position>(r_end);

	bamReader reader(filename, refName, start, end);
	reader.fetch();

    std::vector<aligned_read> const& reads = reader.get_reads();

	field<t_seq_bases> bases(reads.size());
	field<t_epsilon_quality> epsilon(reads.size());
	t_positions pos(reads.size());
	t_lengths len(reads.size());

    std::vector<std::string> name(reads.size());

    for(unsigned int i = 0; i < reads.size(); ++i) {

		aligned_read read = reads[i];

		bases(i) = read.bases;
		epsilon(i) = read.epsilon;
		pos(i) = read.position;
		len(i) = read.length;

        name[i] = read.name;
	}

    rList res;

	res.attach(rObject(bases), "reads");
	res.attach(rObject(epsilon), "quality");
	res.attach(rObject(pos), "positions");
	res.attach(rObject(len), "lengths");
    res.attach(rObject(name), "names");

    return rObject(res);
}

SEXP r_epiG_fetch_reads(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end) {

	try {

		return epiG_fetch_reads(r_filename, r_refName, r_start, r_end);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if (e.what() != NULL) {
			report_error(e.what());
		}

		else {
			report_error("Unknown error");
		}

	} catch (...) {
		report_error("Unknown error");
	}

	return R_NilValue; //Avoid compiler warnings
}

SEXP epiG_fetch_reads_raw(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end) {

	const std::string filename = get_value<std::string>(r_filename);
	const std::string refName = get_value<std::string>(r_refName);
	const t_position start = get_value<t_position>(r_start);
	const t_position end = get_value<t_position>(r_end);

	bamReader reader(filename, refName, start, end);
	reader.fetch_raw();

    std::vector<aligned_read_raw> const& reads = reader.get_reads_raw();

	field<t_seq_bases> bases(reads.size());
	field<t_quality> quality(reads.size());
	t_positions pos(reads.size());
	t_lengths len(reads.size());
    std::vector<std::string> name(reads.size());

    for(unsigned int i = 0; i < reads.size(); ++i) {

		aligned_read_raw read = reads[i];

		bases(i) = read.bases;
		quality(i) = read.quality;
		pos(i) = read.position;
		len(i) = read.length;
        name[i] = read.name;

	}

    rList res;

	res.attach(rObject(bases), "reads");
	res.attach(rObject(quality), "quality");
	res.attach(rObject(pos), "positions");
	res.attach(rObject(len), "lengths");
    res.attach(rObject(name), "names");

    return rObject(res);
}


SEXP r_epiG_fetch_reads_raw(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end) {

	try {

		return epiG_fetch_reads_raw(r_filename, r_refName, r_start, r_end);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if (e.what() != NULL) {
			report_error(e.what());
		}

		else {
			report_error("Unknown error");
		}

	} catch (...) {
		report_error("Unknown error");
	}

	return R_NilValue; //Avoid compiler warnings
}

SEXP epiG_fetch_reads_info(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end) {

    const std::string filename = get_value<std::string>(r_filename);
    const std::string refName = get_value<std::string>(r_refName);
    const t_position start = get_value<t_position>(r_start);
    const t_position end = get_value<t_position>(r_end);

    bamReader reader(filename, refName, start, end);
    reader.fetch_reads_info();
    std::vector<read_info> const& infos = reader.get_infos();

    t_positions pos(infos.size());
    t_lengths len(infos.size());
    std::vector<std::string> name(infos.size());

    for(unsigned int i = 0; i < infos.size(); ++i) {

        read_info info = infos[i];

        pos(i) = info.position;
        len(i) = info.length;
        name[i] = info.name;
    }

    rList res;

    res.attach(rObject(pos), "positions");
    res.attach(rObject(len), "lengths");
    res.attach(rObject(name), "names");

    return rObject(res);
}

SEXP r_epiG_fetch_reads_info(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end) {

    try {

        return epiG_fetch_reads_info(r_filename, r_refName, r_start, r_end);

        //Catch unhandled exceptions

    } catch (std::exception & e) {

        if (e.what() != NULL) {
            report_error(e.what());
        }

        else {
            report_error("Unknown error");
        }

    } catch (...) {
        report_error("Unknown error");
    }

    return R_NilValue; //Avoid compiler warnings
}


SEXP epiG_fetch_header(SEXP r_filename) {

    const std::string filename = get_value<std::string>(r_filename);

    //Open bam file
    samfile_t *file_handle = samopen(filename.c_str(), "rb", 0);

    if (file_handle == 0) {
    	throw std::runtime_error("Fail to open BAM file.\n");
    }

    bam_header_t header = *(file_handle->header);

	std::string txt(header.text, header.l_text);
    std::vector <std::string> targets;
    std::vector <int> target_len;

    int n_targets = header.n_targets;
    for(int i = 0; i < n_targets; i++) {
    	std::string tn(header.target_name[i]);
    	targets.push_back(tn);
    	target_len.push_back(header.target_len[i]);
    }

//    for(int i = 0; i < n_targets; i++) {
//    	std::cout << targets[i] << std::endl;
//    }

    samclose(file_handle);

    rList res;

    res.attach(rObject(targets), "refnames");
    res.attach(rObject(target_len), "lengths");
    res.attach(rObject(txt), "text");

    return rObject(res);
}

SEXP r_epiG_fetch_header(SEXP r_filename) {

    try {

        return epiG_fetch_header(r_filename);

        //Catch unhandled exceptions

    } catch (std::exception & e) {

        if (e.what() != NULL) {
            report_error(e.what());
        }

        else {
            report_error("Unknown error");
        }

    } catch (...) {
        report_error("Unknown error");
    }

    return R_NilValue; //Avoid compiler warnings
}

#endif /* EPIG_FETCH_R_INTERFACE_HPP_ */
