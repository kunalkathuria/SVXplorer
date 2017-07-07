#!/usr/local/bin/perl
 
=head1 NAME
 
    run_cnvnator_on_assembly.pl
 
=head1 SYNOPSIS
 
    run_cnvnator_on_assembly.pl input_fasta input_bam output outputdir path_to_cnvnator windowsize
        where input_fasta is the input fasta file,
              input_bam is the input bam file,
              output is the output file,
              outputdir is the output directory for writing output files,
              path_to_cnvnator is the path to CNVnator,
              windowsize is the window size to use.
 
=head1 DESCRIPTION
 
    This script runs CNVnator on an input assembly, given a fasta file for the
    assembly (<input_fasta>) and a bam file of reads mapped to the assembly
    (<input_bam>).
       Note: before running this script, you need to run:
    setenv ROOTSYS "/nfs/users/nfs_b/bf3/bin/root/"   
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ROOTSYS}/lib
       Note: this job needs to be submitted to the farm with at least 2.5 Gbyte of memory (RAM).
 
=head1 VERSION
  
    Perl script last edited 28-Jan-2013.
 
=head1 CONTACT
 
    alc@sanger.ac.uk (Avril Coghlan)
 
=cut
 
# 
# Perl script run_cnvnator_on_assembly.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 28-Jan-13.
# Last edited 28-Jan-2013.
# SCRIPT SYNOPSIS: run_cnvnator_on_assembly.pl: run CNVnator on an assembly, given the assembly fasta & a bam file of reads mapped to the assembly
#
#------------------------------------------------------------------#
 
# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:
 
use strict;
use warnings;
 
my $num_args               = $#ARGV + 1;
if ($num_args != 6)
{
    print "Usage of run_cnvnator_on_assembly.pl\n\n";
    print "perl run_cnvnator_on_assembly.pl <input_fasta> <input_bam> <output> <outputdir> <path_to_cnvnator> <windowsize>\n";
    print "where <input_fasta> is the input fasta file,\n";
    print "      <input_bam> is the input bam file,\n";
    print "      <output> is the output file,\n";
    print "      <outputdir> is the output directory for writing output files,\n";
    print "      <path_to_cnvnator> is the path to CNVnator,\n";
    print "      <windowsize> is the window size to use\n";
    print "For example, >perl run_cnvnator_on_assembly.pl ref.fa out.sorted.markdup.bam\n";
    print "output /lustre/scratch108/parasites/alc/StrongyloidesCNVnator\n";
    print "/nfs/users/nfs_b/bf3/bin/CNVnator_v0.2.7/src/ 100\n";
    print "Note: before running this script, you need to run:\n";
    print "setenv ROOTSYS \"/nfs/users/nfs_b/bf3/bin/root/\"\n";
    print "setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:\${ROOTSYS}/lib\n";
    exit;
}
 
# FIND THE PATH TO THE INPUT FASTA FILE:                     
 
my $input_fasta            = $ARGV[0];
 
# FIND THE PATH TO THE INPUT BAM FILE:
 
my $input_bam              = $ARGV[1];
 
# FIND THE PATH TO THE OUTPUT FILE:
 
my $output                 = $ARGV[2];
 
# FIND THE DIRECTORY TO USE FOR OUTPUT FILES:      
 
my $outputdir              = $ARGV[3];
 
# FIND THE PATH TO CNVnator:
 
my $path_to_cnvnator       = $ARGV[4];
 
# FIND THE WINDOW SIZE TO USE:
 
my $windowsize             = $ARGV[5];
 
#------------------------------------------------------------------#
 
# TEST SUBROUTINES: 
 
my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
&test_print_error; 
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
&run_main_program($outputdir,$input_fasta,$input_bam,$output,$path_to_cnvnator,$windowsize);
 
print STDERR "FINISHED.\n";
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
sub run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $input_fasta         = $_[1]; # THE INPUT FASTA FILE  
   my $input_bam           = $_[2]; # THE INPUT BAM FILE
   my $output              = $_[3]; # THE OUTPUT FILE 
   my $path_to_cnvnator    = $_[4]; # PATH TO CNVNATOR
   my $windowsize          = $_[5]; # WINDOW SIZE TO USE
   my $errorcode;                   # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg;                    # RETURNED AS 'none' IF THERE IS NO ERROR. 
   my $temp_fasta;                  # A REFORMATTED VERSION OF THE FASTA FILE. 
 
   # REFORMAT THE FASTA FILE SO THAT THERE ARE 60 CHARACTERS PER LINE OF SEQUENCE:
   ($temp_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   system "perl /rhome/cjinfeng/HEG4_cjinfeng/Variations/SV/CNVnator/bin/reformat_fasta.pl $input_fasta $temp_fasta $outputdir";
 
   # RUN CNVNATOR ON ALL SCAFFOLDS IN THE INPUT FASTA FILE:
   ($errorcode,$errormsg)  = &run_cnvnator_on_scaffolds($outputdir,$temp_fasta,$input_bam,$output,$path_to_cnvnator,$windowsize);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   system "rm -f $temp_fasta";
}
 
#------------------------------------------------------------------#
 
# RUN CNVNATOR ON ALL SCAFFOLDS IN THE INPUT FASTA FILE:
 
sub run_cnvnator_on_scaffolds
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES INTO
   my $input_fasta         = $_[1]; # INPUT FASTA FILE FOR ASSEMBLY
   my $input_bam           = $_[2]; # INPUT BAM FILE
   my $output              = $_[3]; # OUTPUT FILE
   my $path_to_cnvnator    = $_[4]; # PATH TO CNVnator
   my $windowsize          = $_[5]; # WINDOW SIZE TO USE
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR 
   my $scaffold;                    # NAME OF A SCAFFOLD 
   my $listfile;                    # FILE WITH THE LIST OF SCAFFOLDS WE WANT THE SEQUENCES FOR 
   my $seqfile;                     # FILE WITH THE SEQUENCE FOR THE SCAFFOLD: 
   my @temp;                        # 
 
   # OPEN THE OUTPUT FILE:
   $output                 = $outputdir."/".$output;
   open(OUTPUT,">$output") || die "ERROR: run_cnvnator_on_scaffolds: cannot open output $output.\n";
   close(OUTPUT);
 
   # GET THE NAMES OF THE SEQUENCES FROM THE FASTA FILE:
   open(TEMP, "grep \">\" $input_fasta | cut -d\">\" -f2 | ");
   while(<TEMP>)
   {
      $scaffold            = $_; 
      chomp $scaffold;
      #next unless ($scaffold eq "chr3");
      print STDERR "================================================================\n";
      print STDERR "Running CNVnator for scaffold $scaffold...\n";
      # MAKE A FILE WITH THE SEQUENCE FOR $scaffold:
      ($listfile,$errorcode,$errormsg) = &make_filename($outputdir);
      print "listfile: $listfile\n";
      if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
      open(LISTFILE,">$listfile") || die "ERROR: run_cnvnator_on_scaffolds: cannot open $listfile\n";
      print LISTFILE "$scaffold\n"; 
      close(LISTFILE);
      # CNVnator REQUIRES THAT YOU HAVE A FASTA FILE FOR EACH SCAFFOLD IN THE CURRENT DIRECTORY:
      $seqfile             = $scaffold.".fa";
      
      print STDERR "Making file $seqfile...\n";
      system "perl /rhome/cjinfeng/HEG4_cjinfeng/Variations/SV/CNVnator/bin/get_seqs_in_list.pl $input_fasta $listfile $outputdir $seqfile";
      $seqfile             = $outputdir."/".$seqfile;
      # WAIT A FEW SECONDS FOR THE SEQUENCE FILE TO BE MADE (OTHERWISE WE SOMETIMES GET ERRORS):
      system "sleep 10";
      if (!(-e $seqfile))
      {
         $errormsg         = "ERROR: run_cnvnator_on_scaffolds: file $seqfile doesn't exist\n";
         $errorcode        = 1; # ERRORCODE=1
         return($errorcode,$errormsg);
      }
      # RUN CNVNATOR FOR THIS SCAFFOLD:
      ($errorcode,$errormsg) = &run_cnvnator_on_scaffold($outputdir,$input_bam,$output,$input_fasta,$scaffold,$path_to_cnvnator,$windowsize,$seqfile);
      if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
      # DELETE THE TEMPORARY FILES:
      system "rm -f $listfile";
      system "rm -f $seqfile";   
   }
   close(TEMP);
 
   return($errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# RUN CNVNATOR ON A SCAFFOLD:
 
sub run_cnvnator_on_scaffold
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES INTO
   my $input_bam           = $_[1]; # INPUT BAM FILE
   my $output              = $_[2]; # OUTPUT FILE
   my $input_fasta         = $_[3]; # INPUT FASTA FILE
   my $scaffold            = $_[4]; # SCAFFOLD NAME
   my $path_to_cnvnator    = $_[5]; # PATH TO CNVnator
   my $windowsize          = $_[6]; # WINDOWSIZE TO USE
   my $seqfile             = $_[7]; # FASTA FILE FOR THE SCAFFOLD
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR 
   my $cmd;                         # COMMAND TO RUN
   my $root;                        # out.root FILE
   my @temp;                        #
 
   # NOTE: SOME SCAFFOLD NAMES HAVE # IN THE SCAFFOLD NAME, SO I HAD TO PUT '$scaffold' IN THESE COMMANDS TO MAKE THEM RUN:
   # RUN -tree
   $root                   = $output.".root";
   @temp                   = split(/\//,$root);
   $root                   = $temp[$#temp];
   print STDERR "Running CNVnator -tree for scaffold $scaffold ...\n";
   $cmd                    = $path_to_cnvnator."/cnvnator -root $root -genome $input_fasta -chrom \'$scaffold\' -tree $input_bam"; 
   print STDERR "Running $cmd...\n";
   system "$cmd";
 
   # RUN -his
   print STDERR "Running CNVnator -his for scaffold $scaffold...\n";
   $cmd                    = $path_to_cnvnator."/cnvnator -root $root -genome $input_fasta -chrom \'$scaffold\' -his $windowsize"; 
   print STDERR "Running $cmd...\n";
   system "$cmd";
 
   # RUN -stat
   print STDERR "Running CNVnator -stat for scaffold $scaffold...\n";
   $cmd                    = $path_to_cnvnator."/cnvnator -root $root -genome $input_fasta -chrom \'$scaffold\' -stat $windowsize";
   print STDERR "Running $cmd...\n";
   system "$cmd";
 
   # RUN -partition 
   print STDERR "Running CNVnator -partition for scaffold $scaffold...\n";
   $cmd                    = $path_to_cnvnator."/cnvnator -root $root -genome $input_fasta -chrom \'$scaffold\' -partition $windowsize"; 
   print STDERR "Running $cmd...\n";
   system "$cmd";
 
   # RUN -call
   print STDERR "Running CNVnator -call for scaffold $scaffold...\n";
   $cmd                    = $path_to_cnvnator."/cnvnator -root $root -genome $input_fasta -chrom \'$scaffold\' -call $windowsize >> $output"; 
   print STDERR "Running $cmd...\n";
   system "$cmd";
 
   # DELETE TEMPORARY FILES:     
   system "rm -f $outputdir/$root";
 
   return($errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# SUBROUTINE TO MAKE A FILE NAME FOR A TEMPORARY FILE:
 
sub make_filename
{
   my $outputdir             = $_[0]; # OUTPUT DIRECTORY TO WRITE TEMPORARY FILE NAME TO
   my $found_name            = 0;     # SAYS WHETHER WE HAVE FOUND A FILE NAME YET
   my $filename              = "none";# NEW FILE NAME TO USE 
   my $errorcode             = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg              = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $poss_filename;                 # POSSIBLE FILE NAME TO USE
   my $random_number;                 # RANDOM NUMBER TO USE IN TEMPORARY FILE NAME
 
   while($found_name == 0)
   {
      $random_number         = rand();
      $poss_filename         = $outputdir."/tmp".$random_number;
      if (!(-e $poss_filename))
      {
         $filename           = $poss_filename;
         $found_name         = 1;
      } 
   } 
   if ($found_name == 0 || $filename eq 'none')
   {
      $errormsg              = "ERROR: make_filename: found_name $found_name filename $filename\n";
      $errorcode             = 6; # ERRORCODE=6 (SHOULDN'T HAPPEN, SO CAN'T TEST)
      return($filename,$errorcode,$errormsg);
   }
 
   return($filename,$errorcode,$errormsg); 
}
 
#------------------------------------------------------------------#
 
# TEST &print_error
 
sub test_print_error
{
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE WAS NO ERROR
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE WAS NO ERROR
 
   ($errormsg,$errorcode)  = &print_error(45,45,1);
   if ($errorcode != 12) { print STDERR "ERROR: test_print_error: failed test1\n"; exit;}
 
   ($errormsg,$errorcode)  = &print_error('My error message','My error message',1);
   if ($errorcode != 11) { print STDERR "ERROR: test_print_error: failed test2\n"; exit;}
 
   ($errormsg,$errorcode)  = &print_error('none',45,1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test3\n"; exit;} 
 
   ($errormsg,$errorcode)  = &print_error('My error message', 0, 1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test4\n"; exit;}
}
 
#------------------------------------------------------------------#
 
# PRINT OUT AN ERROR MESSAGE AND EXIT.
 
sub print_error
{
   my $errormsg            = $_[0]; # THIS SHOULD BE NOT 'none' IF AN ERROR OCCURRED.
   my $errorcode           = $_[1]; # THIS SHOULD NOT BE 0 IF AN ERROR OCCURRED.
   my $called_from_test    = $_[2]; # SAYS WHETHER THIS WAS CALLED FROM test_print_error OR NOT
 
   if ($errorcode =~ /[A-Z]/ || $errorcode =~ /[a-z]/) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 11; $errormsg = "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; # ERRORCODE=11
         return($errormsg,$errorcode);
      }
      else 
      { 
         print STDERR "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; 
         exit;
      }
   }
 
   if (!($errormsg =~ /[A-Z]/ || $errormsg =~ /[a-z]/)) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 12; $errormsg = "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; # ERRORCODE=12
         return($errormsg,$errorcode);
      }
      else
      {
         print STDERR "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; 
         exit;
      }
   }
 
   if    ($errormsg eq 'none' || $errorcode == 0) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 13; $errormsg = "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; # ERRORCODE=13
         return($errormsg,$errorcode);
      }
      else 
      {
         print STDERR "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; 
         exit;
      }
   }
   else                                           
   { 
      print STDERR "$errormsg"; 
      exit;                                                      
   } 
 
   return($errormsg,$errorcode);
}
 
#------------------------------------------------------------------#

