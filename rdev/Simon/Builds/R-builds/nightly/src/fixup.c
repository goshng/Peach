#include <unistd.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  setuid(0);
  system("/usr/sbin/chown -R root:admin /Library/Frameworks/RNAseq.framework");
  system("/bin/chmod -R g+w /Library/Frameworks/RNAseq.framework");
  return 0;
}
