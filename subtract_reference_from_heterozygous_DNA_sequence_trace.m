function subtract_reference_from_heterozygous_DNA_sequence_trace()
%SUBTRACT_REFERENCE_FROM_HETEROZYGOUS_DNA_SEQUENCE_TRACE  Subtract a
%reference sequence from a heterozygous sequence trace in order to get the
%sequence of the homologous chromosome
% SUBTRACT_REFERENCE_FROM_HETEROZYGOUS_DNA_SEQUENCE_TRACE
% 
% For sequence traces of heterozygous parts of the genome, some positions
% will exhibit "mixed peaks" in which two nucleotides are dominant, since
% the two chromosomes differ in sequence at that position. It's easy enough
% to detect these locations manually, by examining the sequence trace.
% However, this becomes challenging if one of the chromosomes contains an
% insertion/deletion (indel), in which case all positions beyond the
% location of the indel will exhibit mixed peaks, making it very difficult
% to subtract the reference sequence by manual examination.
%
%
% USAGE
%
% To use this function, first convert sequence files from to .fcs. If they
% are currently .ab1 and the user is using OS X, they can use the following
% instructions (modified from this website:
% https://www.biostars.org/p/622/):
% 1. Open terminal
% 2. cd to progs
% 3. Run convert_trace command; for example:
% ./convert_trace -out_format scf < salr_EagI_PCR_rev-PAA-04.ab1 > salr_EagI_PCR_rev-PAA-04.fcs
% 4. move output file from progs directory to desired directory
%
% This function doesn't have any input or output arguments. Simply type the
% reference and sample filenames (with the full path, if the files are in a
% different directory as this source code) below in the first lines, and 
% run the program. The program prints the amount that the sample
% sequence needed to be shifted in order to align it with the reference
% sequence, and then it prints the sample sequence after the reference
% sequence has been subtracted from it. To detect the presence of indels,
% one can simply align this output to the reference sequence (using tools
% such as BLAST).
%
% In cases where the subtracted sequence is ambiguous, IUPAC conventions
% are used to indicate the level of ambiguity.
%
% The program uses basepairs 61-80 to align the two sequences. Thus, the 
% two sequences should be at least ~100 basepairs long.
%
%
% IMPLEMENTATION
%
% First, the two sequences are aligned using basepairs 61-80 of the sample
% sequences (we didn't want to use an earlier part of the sequence since
% sequencing tends to be less reliable at the start of the target). Then,
% the reference sequence is subtracted from the sample sequence.
%
% ******
% Created by Abed Alnaif, abed.alnaif@gmail.com
% Tested in Matlab R2012b
% ******

%% ***** USER INPUTS *****
% Specify the locations of sequence files
% specify filename containing trace for reference sequence
reference_filename = 'salr_EagI_PCR_rev-PAA-04.fcs';

% specify filename containing trace for sample sequence
sample_filename = 'U290_salr_EagI-EagI_for.fcs';

% ***** END OF USER INPUTS *****

% print 'sample_filename' to screen
disp(['reference: ' reference_filename]);
disp(['sample: ' sample_filename]);

% read in data
[reference.sample,reference.probability] = scfread(reference_filename);
[sample.sample,sample.probability] = scfread(sample_filename);

% align the sequences using nucleotides 61-80 of sample
equiv_index = strfind(reference.probability.base',sample.probability.base(61:80)');

% determine how much to shift reference sequence
shift = equiv_index - 61;

disp(['shift to align reference and sample traces: ' num2str(shift)])

% do not perform any calculations on first 60 nucleotides
sequence_with_reference_subtracted = sample.probability.base(1:60)';

% specify ambiguous nucleotides by 'N' instead of '-', so as not to
% conflict with BLAST, which uses '-' to indicate deletions
sequence_with_reference_subtracted = strrep(sequence_with_reference_subtracted,'-','N'); 

nucleotides = {'A','C','G','T'};

% loop through all nucleotides of sample and subtract reference
nNucleotidesSample = length(sample.probability.peak_index);
nNucleotidesRef = length(reference.probability.peak_index);
for sample_index=61:nNucleotidesSample % loop through all peaks
    % debug lines
%     if sample_index == 189
%         1;
%     end
    
    % determine location of peak
    peak_location_sample = sample.probability.peak_index(sample_index);
    
    % find nucleotide with max peak
    sample_peak_val = nan(4,1);
    sample_peak_val(1) = sample.sample.A(peak_location_sample);
    sample_peak_val(2) = sample.sample.C(peak_location_sample);
    sample_peak_val(3) = sample.sample.G(peak_location_sample);
    sample_peak_val(4) = sample.sample.T(peak_location_sample);
    
    % determine size of max peak and nucleotide corresponding to max peak
    [max_val_sample,max_index_sample] = max(sample_peak_val);
    
    output_nucleotide_indices = NaN;
    
    % define threshold as 100, otherwise peaks are poor quality
    if max_val_sample < 100
        output_nucleotide_indices = NaN;
    else

        % determine how many nucleotides exhibit peak at this location,
        % defining threshold as max/5
        sample_nucleotide_indices_in_peak = find(sample_peak_val > max_val_sample/5);
        sample_num_peaks = length(sample_nucleotide_indices_in_peak);

        if sample_num_peaks == 1
            % there's only one nucleotide in this peak; no need to subtract
            % reference
            output_nucleotide_indices = sample_nucleotide_indices_in_peak;
        elseif sample_num_peaks > 2
            % if more than two nucleotides in peak, probably the data is
            % poor quality
            output_nucleotide_indices = NaN;
        else
            % there are two peaks; subtract reference
            reference_index = sample_index + shift; % determine equivalent index in reference trace
            
            if reference_index < 1 || reference_index > nNucleotidesRef
                % no data for reference since we are outside the bounds of
                % sequence; thus, I cannot subtract reference
                output_nucleotide_indices = sample_nucleotide_indices_in_peak;
            else
                peak_location_reference = reference.probability.peak_index(reference_index);
                ref_peak_val = nan(4,1);
                ref_peak_val(1) = reference.sample.A(peak_location_reference);
                ref_peak_val(2) = reference.sample.C(peak_location_reference);
                ref_peak_val(3) = reference.sample.G(peak_location_reference);
                ref_peak_val(4) = reference.sample.T(peak_location_reference);
                [max_val_ref,max_index_ref] = max(ref_peak_val);

                if max_val_ref < 100
                    % reference data is poor quality; cannot subtract
                    % reference; take peak to be either of the two peaks in
                    % sample trace
                    output_nucleotide_indices = sample_nucleotide_indices_in_peak;
                else
                    ref_nucleotide_indices_in_peak = find(ref_peak_val > max_val_ref/5);
                    ref_num_peaks = length(ref_nucleotide_indices_in_peak);
                    if ref_num_peaks == 1
                        % subtract nucleotide corresponding to reference
                        output_nucleotide_indices = sample_nucleotide_indices_in_peak(sample_nucleotide_indices_in_peak~=ref_nucleotide_indices_in_peak);
                    elseif ref_num_peaks > 2
                        % reference exhibits mixed peaks; cannot subtract reference
                        output_nucleotide_indices = sample_nucleotide_indices_in_peak;
                    end
                end
            end
        end
    end
    
    if (length(output_nucleotide_indices)==1)
        % there is only one index
        if isnan(output_nucleotide_indices)
            output_nucleotide = 'N';
        else
            % a single nucleotide was successfully identified
            output_nucleotide = nucleotides{output_nucleotide_indices};
        end
    else
        % there are two indices (two possible nucleotides); use IUPAC code
        if output_nucleotide_indices == [1;2]
            % A or C
            output_nucleotide = 'M';
        elseif output_nucleotide_indices == [1;3]
            % A or G
            output_nucleotide = 'R';
        elseif output_nucleotide_indices == [1;4]
            % A or T
            output_nucleotide = 'W';
        elseif output_nucleotide_indices == [2;3]
            % C or G
            output_nucleotide = 'S';
        elseif output_nucleotide_indices == [2;4]
            % C or T
            output_nucleotide = 'Y';
        elseif output_nucleotide_indices == [3;4]
            % G or T
            output_nucleotide = 'K';
        end
    end
    
    sequence_with_reference_subtracted = [sequence_with_reference_subtracted output_nucleotide];
end

disp(sequence_with_reference_subtracted)

end