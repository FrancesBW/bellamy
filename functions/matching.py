#! /usr/bin/env python

def reference_pre_filter(target_catalog, reference_catalog, target_freq=None, reference_freq=None):
        """
        Read the fits files for the two catalogues, the target and the reference catalogue and returns a subset of the reference catalog. The reference catalogue is likely bigger than the target catalogue and may contain multi-frequency data, this function reduces the amount of data searched later on by throwing away all sources in the reference catalogue outside the limits of the target catalogue (with a crude buffer of 5*max_psf of the target catalogue) and in futrue versions will also throw away data from frequencies not covered in the target catalogue. Future versions should also determine the resolution of the reference catalog to determine the edge parameter.
        
        :param target_catalog: filename for the target catalogue
        :param reference_catalog: filename for the reference catalogue

        :param target_freq: frequencies covered by the target catalogue formatted as list or np.array (currently not utilised,optional)
        :param reference_freq: frequencies covered by the reference catalogue formatted as list or np.array (currently not utilised,optional)

        :return: updated_reference_catalog
        """

        #determining the maximum and minimum ra and dec of the target catalog
        max_ra=np.max(target_catalog['ra'])
        min_ra=np.min(target_catalog['ra'])
        max_dec=np.max(target_catalog['dec'])
        min_dec=np.min(target_catalog['dec'])
        #determining the maximum source half width in the target catalog
        max_a=np.max(target_catalog['a'])
        max_b=np.max(target_catalog['b'])
        #determining the maximum beam half width in the target catalog
        max_beam_a=np.max(target_catalog['psf_a'])
        max_beam_b=np.max(target_catalog['psf_b'])
        #adding these both together to give generous estimate of offset allowance for the target catalogue from the reference catalogue
        edge=max(max_a,max_b)/3600+max(max_beam_a,max_beam_b)/3600

        #setting the maximum and minimum ra and dec to be pulled from the reference catalog with 5* the offset allowance previously calculated
        filter_max_ra=max_ra+5*edge
        filter_min_ra=min_ra-5*edge
        filter_max_dec=max_dec+5*edge
        filter_min_dec=min_dec-5*edge
        #setting the filter with the conditions
        reference_filter=np.where((reference_catalog['ra']<filter_max_ra) & (reference_catalog['ra']>filter_min_ra) & (reference_catalog['dec']<filter_max_dec) & (reference_catalog['dec']>filter_min_dec))
        return reference_catalog[reference_filter]



def model_offsets_and_update_positions(cross_matched_catalogue,target_catalogue,run_num, ref_ra_column='ref_ra', ref_dec_column='ref_dec', tar_ra_column='tar_ra', tar_dec_column='tar_dec'):
        """
        Read in the cross matched catalogue and model the measured offsets.
        
        :param cross_matched_catalogue: cross matched catalogue in Table format, output from 'cross_matching' function will provide correct format to use default optional parameters
        :param target_catalogue: Table format of target catalogue with image positions (unupdated by previous models)
        :param run_num: Iteration identifier to allow for models and offsets at each step of the process to be saved as separate image files

        :param ref_ra_column: the column value for the source position ra in the reference catalogue (optional)
        :param ref_dec_column: the column value for the source position dec in the reference catalogue (optional)
        :param tar_ra_column: the column value for the source position ra in the target catalogue (optional)
        :param tar_dec_column: the column value for the source position dec in the target catalogue (optional)

        :return:  cross_matched_table containing position information of matched sources only, reference catalogue with cross matched sources removed, modelled position target catalogue with cross matched sources removed, original position target catalogue with cross matched sources removed
        """

        #calculating the offsets from the cross matched catalogue
        d_ra=cross_matched_catalogue[tar_ra_column]-cross_matched_catalogue[ref_ra_column]
	d_dec=cross_matched_catalogue[tar_dec_column]-cross_matched_catalogue[ref_dec_column]
        
        #creating a model of the measured offsets as a function of the position in the sky using rbf interpolation
	model_d_ra=interpolate.Rbf(cross_matched_catalogue[tar_ra_column], cross_matched_catalogue[tar_dec_column], d_ra, function='linear',smooth=0.032777778)
	model_d_dec=interpolate.Rbf(cross_matched_catalogue[tar_ra_column], cross_matched_catalogue[tar_dec_column], d_dec, function='linear',smooth=0.032777778)
        
        #need to create a copy of the target catalogue to avoid overwriting the original positions at this point
        target_catalogue_copy=copy(target_catalogue)
        
        #updating the target catalogue to try to bring it in line with the GLEAM catalogue
        target_catalogue_copy['ra']=target_catalogue_copy['ra']-model_d_ra(target_catalogue_copy['ra'], target_catalogue_copy['dec'])
        target_catalogue_copy['dec']=target_catalogue_copy['dec']-model_d_dec(target_catalogue_copy['ra'], target_catalogue_copy['dec'])
        #print(Table([target_catalogue['ra'],target_catalogue['dec'],target_catalogue['uuid']],names=('tar_ra','tar_dec','uuid')))

        modelled_offset_ra=d_ra-model_d_ra(cross_matched_catalogue[tar_ra_column],cross_matched_catalogue[tar_dec_column])
        modelled_offset_dec=d_dec-model_d_dec(cross_matched_catalogue[tar_ra_column], cross_matched_catalogue[tar_dec_column])
        #print(d_ra)
        #print(model_d_ra(cross_matched_catalogue[tar_ra_column],cross_matched_catalogue[tar_dec_column]))        
        if g_plot==True:
                #making quiver plot to show measured offsets
                fig = plt.figure(figsize=(12, 12))
                gs = gridspec.GridSpec(100,100)
                gs.update(hspace=0,wspace=0)
                ax = fig.add_subplot(gs[0:100,0:100])
                angles = np.degrees(np.arctan2(d_dec, d_ra))
                cax = ax.quiver(cross_matched_catalogue[tar_ra_column], cross_matched_catalogue[tar_dec_column], d_ra, d_dec,angles,cmap=plt.cm.get_cmap('rainbow'))
                ax.set_xlabel("Distance from pointing centre / degrees")
                ax.set_ylabel("Distance from pointing centre / degrees")
                ax.set_title("Source position offsets / arcsec")
                plt.savefig("Measured_offsets" +g_file_suffix+"_run_"+str(run_num)+".png")
                plt.close()

                #making quiver plot to show uniformly sampled modelled offsets
                fig = plt.figure(figsize=(12,12))
                gs = gridspec.GridSpec(100,100)
                gs.update(hspace=0,wspace=0)
                ax = fig.add_subplot(gs[0:100,0:100])
                #gx=np.linspace(70.0,99.0,1000)
                #gy=np.linspace(-17.0,11.0,1000)
                gx, gy = np.mgrid[70.:99.:(29.0)/50.,-17.:11.:(28.0)/50.]
                angles = np.degrees(np.arctan2(model_d_dec(gx,gy), model_d_ra(gx,gy)))
                cax = ax.quiver(gx, gy, model_d_ra(gx, gy), model_d_dec(gx,gy),angles,cmap=plt.cm.get_cmap('rainbow'))
                ax.set_xlabel("Distance from pointing centre / degrees")
                ax.set_ylabel("Distance from pointing centre / degrees")
                ax.set_title("Source position offsets / arcsec")
                plt.savefig("Modelled_offsets" +g_file_suffix+"_run_"+str(run_num)+".png")
                plt.close()
        
        return target_catalogue_copy

def flux_prob(ref_candidates, tar_candidates):
        """
        To calculate the gaussian probability that a source matches in flux. Altered positions must be passed to this function, not original ones.
        
        :param ref_candidates: the reference sources that may match the target candidates
        :param tar_candidates: the target source being matched

        :return:  an array that represents the relative probability of the target source matching each reference source
        """
        ultimate_error_allowance=np.sqrt(ref_candidates['local_rms']**2+tar_candidates['local_rms']**2+tar_candidates['err_peak_flux']**2+ref_candidates['err_peak_flux']**2)
        true_error=ref_candidates['peak_flux']-tar_candidates['peak_flux']
        probs=np.exp(-((true_error)**2)/(2*ultimate_error_allowance**2))
        return probs

def position_prob(ref_candidates, tar_candidates):
        """
        To calculate the gaussian probability that a source matches in position. Altered positions must be passed to this function, not original ones.
        
        :param ref_candidates: the reference sources that may match the target candidates
        :param tar_candidates: the target source being matched

        :return:  an array that represents the relative probability of the target source matching each reference source
        """
        
        ref_resolution=ref_candidates['psf_a']/3600
        ref_position_error=np.maximum(ref_candidates['err_ra'],ref_candidates['err_dec'])
        tar_resolution=tar_candidates['psf_a']/3600
        tar_position_error=np.maximum(tar_candidates['err_ra'],tar_candidates['err_dec'])
        ultimate_error_allowance=np.sqrt(ref_resolution**2+ref_position_error**2+tar_resolution**2+tar_position_error**2)
        true_error=np.sqrt((ref_candidates['ra']-tar_candidates['ra'])**2+(ref_candidates['dec']-tar_candidates['dec'])**2)
        probs=np.exp(-((true_error)**2)/(2*ultimate_error_allowance**2))
        return probs



def prob_comb(ref_candidates, tar_entry, confidence_percentile, flux_on ,final_run):
        """
        To calculate the normalised probabilty of a match on position and flux (optional) and return a confirmed match above some probability threshold.
        
        :param ref_candidates: the reference sources that may match the target candidates
        :param tar_entry: the target source being matched
        :param confidence_percentile: float between 0 and 1 to indicate the required confidence interval for a match
        :param flux_on: condition for inclusion of flux in determining of a match. If True, flux probability will be considered in calculation.
        :param final_run: Flags that program is on the final iteration and ignores the confidence percentile, instead returning a match for every target source

        :return: If a match has been made, returns array with 5 elements, the target source uuid, the matched reference source name, raw probability of a match, normalised probability of a match (NaN if there was only one candidate) and number of candidates. If a match is not made, then False is returned. If final run condition is True, a match is made by default.
        """
        
        pos_prob=position_prob(ref_candidates, tar_entry)
        if flux_on==True:
                flux_probs=flux_prob(ref_candidates, tar_entry)
                combined_probs=pos_prob*flux_probs
        else:
                combined_probs=pos_prob
        normalisation_factor=np.sum(combined_probs)
        final_probs=combined_probs/normalisation_factor
        rounded_probs=final_probs.round(decimals=2)
        negligible_filter=np.where(rounded_probs!=0.)
        combined_probs=combined_probs[negligible_filter]
        names=ref_candidates['Name'][negligible_filter]
        #print(Table(zip(pos_prob,flux_probs,combined_probs,final_probs),names=names))
        table=zip(names,final_probs[negligible_filter])
        probability_table=[tar_entry['uuid'],table]

        single_candidate_confidence=g_single_match_percentile
        
        if final_run==True:
                confidence_percentile=0.0
                single_candidate_confidence=0.0
                
        if len(ref_candidates[negligible_filter])>1:
                most_likely=np.argmax(np.array(probability_table[1])[:,1])
                if float(np.array(probability_table[1])[most_likely,1])>=confidence_percentile:
                        return [probability_table[0],np.array(probability_table[1])[most_likely,0],combined_probs[most_likely],float(np.array(probability_table[1])[most_likely,1]),len(ref_candidates[negligible_filter])]
                else:
                        return False
        else:
                if combined_probs>=single_candidate_confidence:
                        return [probability_table[0],np.array(probability_table[1])[0,0],combined_probs[0],float('NaN'),1]
                else:
                        return False
                        

def cross_matching(ref_catalogue, pre_snr_tar_catalogue, original_dist_tar_catalogue, confidence_percentile, snr_restriction=False, flux_match=True,final_run=False):
        """
        Take the updated reference and target catalogues to perform a cross match within a resolution radius. Allows user to define a restriction on the SNR threshold of source, how tightly the fluxes should match and a normalisation factor for fluxes between the two catalogues.
        
        :param ref_catalogue: the reference catalogue with only un-cross-matched sources remaining in it
        :param pre_snr_tar_catalogue: the target catalogue with only un-cross-matched source remaining in it, which hasn't been filtered by SNR yet.
        :param original_dist_tar_catalogue: same sources as pre_snr_tar_catalogue but with original positions of all target sources
        :param confidence_percentile: float between 0 and 1 to indicate the required confidence interval for a match

        :param snr_restriction: an integer or float that defines some lower limit on SNR. Only target sources above this SNR will be cross matched (optional, default is False)
        :param flux_match: If this condition is True matching probability will be altered by how close the sources match in flux (optional, default is True)
        :param final_run: argument required for comb_prob to override the set confidence percentile

        :return: target catalogue with modelled offsets applied
        """
        
        if final_run==True:
                print('Returning most likely match to remaining '+str(len(pre_snr_tar_catalogue))+' unmatched sources')
        else:
                print('Initialising match of '+str(len(pre_snr_tar_catalogue))+' target sources')
        #applying snr filter if it is specified 
        if snr_restriction!=False:
                snr_filter=np.where((pre_snr_tar_catalogue['peak_flux']/pre_snr_tar_catalogue['local_rms'])>=snr_restriction)
                tar_catalogue=pre_snr_tar_catalogue[snr_filter]
        else:
                tar_catalogue=pre_snr_tar_catalogue

        #we're going to cast a wide net of 10' to find multiple matches per target source
        limiting_res=600.0/3600.0*u.degree

        #Using SkyCoord package to prepare the ref and tar catalogues for the cross match
        ref_cat=SkyCoord(ref_catalogue['ra'], ref_catalogue['dec'], unit=u.degree, frame='icrs')
        tar_cat=SkyCoord(tar_catalogue['ra'], tar_catalogue['dec'], unit=u.degree, frame='icrs')

        ref_cat_uuid=[]
        tar_cat_uuid=[]
        cross_matched_table=Table(names=('tar_ra', 'tar_dec','tar_a','tar_b','tar_pa', 'tar_flux','ref_ra', 'ref_dec','ref_a','ref_b','ref_pa','ref_flux','tar_uuid','ref_name','raw_prob','norm_prob','num_of_candidates'),dtype=('f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','U36','U14','f8','f8','i4'))
        gross_matched_idx_ref, gross_matched_idx_tar, gross_matched_sep, dum=tar_cat.search_around_sky(ref_cat,limiting_res)
        if snr_restriction!=False:
            tar_cat_matched_within_resolution=original_dist_tar_catalogue[snr_filter][gross_matched_idx_tar]
        else:
            tar_cat_matched_within_resolution=original_dist_tar_catalogue[gross_matched_idx_tar]

        #this is becuse the probability determination uses the adjusted positions
        tar_cat_for_prob=tar_catalogue[gross_matched_idx_tar]
        
        ref_cat_matched_within_resolution=ref_catalogue[gross_matched_idx_ref]
        
        for item in tar_cat_for_prob:
                if np.any(np.isin(tar_cat_uuid,item['uuid']))==False:
                        reference_matches=ref_cat_matched_within_resolution[np.where(tar_cat_for_prob['uuid']==item['uuid'])]
                        match=prob_comb(reference_matches, item, confidence_percentile, flux_match,final_run)
                        if match!=False:
                                tar_uuid=match[0]
                                ref_uuid=match[1]
                                target_entry=original_dist_tar_catalogue[np.where(original_dist_tar_catalogue['uuid']==tar_uuid)]
                                reference_entry_idx=np.where(ref_catalogue['Name']==ref_uuid)
                                cross_matched_table.add_row((target_entry['ra'],target_entry['dec'],target_entry['a'],target_entry['b'],target_entry['pa'],target_entry['peak_flux'],ref_catalogue['RAJ2000'][reference_entry_idx],ref_catalogue['DEJ2000'][reference_entry_idx],ref_catalogue['a_181'][reference_entry_idx],ref_catalogue['b_181'][reference_entry_idx],ref_catalogue['pa_181'][reference_entry_idx],np.mean((ref_catalogue['peak_flux_181'],ref_catalogue['peak_flux_189']),axis=0)[reference_entry_idx],target_entry['uuid'],ref_catalogue['Name'][reference_entry_idx],match[2],match[3],match[4]))
                                ref_cat_uuid.append(ref_uuid)
                                tar_cat_uuid.append(tar_uuid)
                                
        ref_cat_uuid=np.array(ref_cat_uuid)
        tar_cat_uuid=np.array(tar_cat_uuid)

                        
        #remove the correctly cross matched sources from the reference catalogue 
        new_ref_index_list=np.delete(np.arange(0,ref_catalogue.shape[0]),np.arange(0,ref_catalogue.shape[0])[np.isin(ref_catalogue['uuid'],ref_cat_uuid)])
        new_ref_catalogue=ref_catalogue[new_ref_index_list]

        #remove the correctly cross matched sources from the modelled position target catalogue and original position target catalogue
        new_tar_index_list=np.delete(np.arange(0,pre_snr_tar_catalogue.shape[0]),np.arange(0,pre_snr_tar_catalogue.shape[0])[np.isin(pre_snr_tar_catalogue['uuid'],tar_cat_uuid)])
        new_tar_catalogue=pre_snr_tar_catalogue[new_tar_index_list]
        new_original_dist_tar_catalogue=original_dist_tar_catalogue[new_tar_index_list]
        
        return cross_matched_table, new_ref_catalogue, new_tar_catalogue, new_original_dist_tar_catalogue

def flux_model(raw_tar_catalogue,filtered_ref_catalogue):
        """
        Take the raw target and reference catalogue, do a quick and dirty cross match of the two with only high snr (10) target sources. Construct a 2D polynomial model of the flux correction factor across the field of sources.
        
        :param raw_tar_catalogue: the target catalogue
        :param filtered_ref_catalogue: the reference catalogue filtered for sources within the RA and DEC range of the target catalogue

        :return: flux correction factor model to apply to target catalogue
        """

        #filter for snr>10 sources
        snr_filter=np.where((raw_tar_catalogue['peak_flux']/raw_tar_catalogue['local_rms'])>=10)
        tar_catalogue=raw_tar_catalogue[snr_filter]

        #we want quite a tight match
        limiting_res=118.0/3600.0*u.degree
        ref_cat=SkyCoord(filtered_ref_catalogue['ra'], filtered_ref_catalogue['dec'], unit=u.degree, frame='icrs')
        tar_cat=SkyCoord(tar_catalogue['ra'], tar_catalogue['dec'], unit=u.degree, frame='icrs')
        gross_matched_idx, gross_matched_sep, dum=tar_cat.match_to_catalog_sky(ref_cat)
        tar_cat_matched_within_resolution=tar_catalogue[np.where(gross_matched_sep<=limiting_res)]
        ref_cat_matched_within_resolution=filtered_ref_catalogue[gross_matched_idx][np.where(gross_matched_sep<=limiting_res)]

        #get our fluxes and ratios for the target and reference catalogues
        reference_flux=ref_cat_matched_within_resolution['peak_flux']
        flux_ratio=tar_cat_matched_within_resolution['peak_flux']/reference_flux

        #use a polynomial 2D model for the field of sources
        resulting_model=interpolate.LSQBivariateSpline(tar_cat_matched_within_resolution['ra'],tar_cat_matched_within_resolution['dec'],flux_ratio,[min(tar_cat_matched_within_resolution['ra']),max(tar_cat_matched_within_resolution['ra'])],[min(tar_cat_matched_within_resolution['dec']),max(tar_cat_matched_within_resolution['dec'])],kx=g_flux_model_deg,ky=g_flux_model_deg)

        if g_plot==True:
                #plot resulting model for user
                x=tar_cat_matched_within_resolution['ra']
                y=tar_cat_matched_within_resolution['dec']
                xi = np.linspace(min(x),max(x),100)
                yi = np.linspace(min(y),max(y),100)
                zi = interpolate.griddata((x, y), flux_ratio, (xi[None,:], yi[:,None]), method='cubic')
                x = np.linspace(min(tar_catalogue['ra']),max(tar_catalogue['ra']),100)
                y = np.linspace(min(tar_catalogue['dec']),max(tar_catalogue['dec']),100)
                X, Y = np.meshgrid(x,y)
                solution=resulting_model.ev(X,Y)
                solution_masked=np.ma.masked_array(solution,mask=np.isnan(zi))
                plot=plt.imshow(solution_masked)
                plt.colorbar(plot)
                plt.savefig('Flux_correction_model.png')
                plt.close()
        return resulting_model

def run():
        """
        Main algorithm to run the iterations of cross matching. Takes no arguments, uses global variables only.

        :return: None
        """
        
        target_sources_num=len(raw_target_table)
        #pre filter the reference catalogue for increased efficiency
        filtered_GLEAM=reference_pre_filter(raw_target_table, raw_reference_table)
        print('The filtered reference catalogue has '+str(len(filtered_GLEAM))+' entries')


        model_flux=options.model_flux

        if model_flux==True:
                model=flux_model(raw_target_table,filtered_GLEAM)
                raw_target_table['peak_flux']=raw_target_table['peak_flux']/model.ev(raw_target_table['ra'],raw_target_table['dec'])
                raw_target_table['err_peak_flux']=raw_target_table['err_peak_flux']/model.ev(raw_target_table['ra'],raw_target_table['dec'])
                del(model)

        #run initial cross match
        print('Run 1')
        cross_match_table,updated_ref_cat,updated_tar_cat,updated_tar_cat_orig_dist=cross_matching(filtered_GLEAM, raw_target_table, raw_target_table, g_multiple_match_percentile, snr_restriction=snr_restrict,flux_match=g_flux_match)

        #modelling will fail if we have less than 2 points, so raise an error before scipy does
        if len(cross_match_table)<2:
                log.error("Cannot create model with 1 match or less. Consider relaxing matching percentile restraints or double check column units are correct")
                sys.exit(1)
                

        #set iteration count
        count=1
        improved=True

        #loop which runs through the process of modelling and cross matching until no improvement is found
        while improved==True:
                print('Run '+str(count+1))
                print('Number of cross matches so far: '+str(len(cross_match_table)))
                start_of_run_cross_match_num=len(cross_match_table)
                adjusted_tar_cat=model_offsets_and_update_positions(cross_match_table,updated_tar_cat_orig_dist,count)
                additional_cross_matches,updated_ref_cat,updated_tar_cat,updated_tar_cat_orig_dist=cross_matching(updated_ref_cat,adjusted_tar_cat, updated_tar_cat_orig_dist, multiple_match_percentile,flux_match=g_flux_match)
                #add the new cross matches to the total table
                cross_match_table=vstack([cross_match_table,additional_cross_matches])
                end_of_run_cross_match_num=len(cross_match_table)
                if end_of_run_cross_match_num==start_of_run_cross_match_num:
                        less_certain_cross_matches,updated_ref_cat,updated_tar_cat,updated_tar_cat_orig_dist=cross_matching(updated_ref_cat,adjusted_tar_cat, updated_tar_cat_orig_dist, g_multiple_match_percentile,flux_match=g_flux_match,final_run=True)
                        #add the new cross matches to the total table
                        cross_match_table=vstack([cross_match_table,less_certain_cross_matches])
                        improved=False
                count+=1
        #determine number of target sources that still are not cross matched
        leftover_tar_sources=len(updated_tar_cat)

        print('Matched '+str(int(target_sources_num-leftover_tar_sources))+ ' out of '+str(int(target_sources_num))+' sources in target catalogue')

        Table(updated_ref_cat).write('leftover_reference_catalogue'+g_file_suffix+'.'+output_format, format=g_output_format,overwrite=True)
        Table(updated_tar_cat_orig_dist).write('leftover_target_catalogue'+g_file_suffix+'.'+output_format, format=g_output_format,overwrite=True)
        Table(updated_tar_cat).write('leftover_unwarped_target_catalogue'+g_file_suffix+'.'+output_format, format=g_output_format,overwrite=True)
        cross_match_table.write('cross_matched_table'+g_file_suffix+'.'+output_format, format=g_output_format,overwrite=True)