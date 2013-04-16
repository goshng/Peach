my @array=(1,2,2,2,1,0);
print "Original length: ";
print scalar (@array);
print "\n";
@array=grep { $_ > 1 } @array;
print join(" xxx ", @array);
print "\n";
print "Length: ";
print scalar(@array);
print "\n";
# $#new_array+1
