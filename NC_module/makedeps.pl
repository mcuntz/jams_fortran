#!/usr/bin/perl
#
# Usage: makedeps.pl [OBJPATH] [SRCPATH]
#
# Generate Makefile object dependencies
# To run makedeps.pl, it will be necessary to modify the first line of 
# this script to point to the actual location of Perl on your system.
#
# Matthias Cuntz, Sep. 2010
# Modified from Uwe Schulzweida's createMakefiles.pl for Echam4
#
$PROG=`basename "$0"`;
#
# Create Makefile
#
open(MAKEFILE, "> make.deps");
print "create new make.deps\n";
print MAKEFILE "# -*- Makefile -*- \n";
print MAKEFILE "# Generated automatically by $PROG \n";
# Dependency listings
$objpath = "";
$srcpath = ".";
if (@ARGV >= 1) {
  $objpath = "$ARGV[0]/";
}
if (@ARGV >= 2) {
  $srcpath = "$ARGV[1]/";
}
&MakeDependsf90($ARGV[2]);
#
exit;
#
#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   chdir($srcpath);
   foreach $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f90$/.o/;
         }
      }
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.f90>) {
      ($objfile = $file) =~ s/\.f90$/.o/;
      open(FILE, $file);
      while (<FILE>) {
         /^\s*#*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
      }

      if (defined @incs) {
          $jn = $#incs;
          while ($jn >= 0) {
             if (-e "../include/$incs[$jn]") {
                $incs[$jn] = "\$(INCLUDE)/$incs[$jn]";
             }
             else {
                print "  $file: include file $incs[$jn] not in \$(INCLUDE)\n";
                pop(@incs);
             }
             $jn--;
          }
      }
      if (defined @incs && $#incs < 0) {undef @incs;}

      if (defined @modules) {
          $jn = $#modules;
          while ($jn >= 0) {
             if ("$objfile" eq "$filename{$modules[$jn]}") {
                #print "  $objfile: $modules[$jn]\n";
                pop(@modules);
             }
             $jn--;
          }
      }
      if (defined @modules && $#modules < 0) {undef @modules;}

      if (defined @incs || defined @modules) {
         print MAKEFILE "$objpath$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
	   if ("$filename{$module}" ne "") {
	     push(@dependencies, "$objpath$filename{$module}");
	   }
	 }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         }
      else {
         print MAKEFILE "$objpath$objfile:\n";
         }
      }
   }
