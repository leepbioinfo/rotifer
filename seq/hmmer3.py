#!/usr/bin/env python3

__version__ = "0.01"

### Import rotifer package
import sys
import os
sys.path.insert(0, os.path.join(os.getcwd(), '../..'))

### Import core cli
import rotifer.core.cli as corecli

class hmmerparser:
    def __init__(self,input_file):
        # self.input_file = input_file
        if isinstance(input_file, io.IOBase):
            self.input_file = input_file
        else:
            self.input_file = open(input_file)
    def parse(self):
        pass
'''
Domain annotation for each model (and alignments):
>> Ammonium_transp  Ammonium Transporter Family
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  287.1  24.8   6.9e-89   1.6e-85       1     395 [.      10     381 ..      10     384 .. 0.95

  Alignments for each domain:
  == domain 1  score: 287.1 bits;  conditional E-value: 6.9e-89
                      HHHHHHHHHHHHHHHHHHHHHHHHSSGGGHHHHHHHHHHHHHHHHHHHHHTHHHHHHSSCBSSSSSCE-TTTTSCCGCCTTCCHHHHHHHHHHH CS
  Ammonium_transp   1 aflllsaalvifmqagfalleaglvrsknvlnilvknlldlavvvllyvlfGyslafgkskgvsgfignlglsaagvqdetlldglfflfqlaf 94 
                      ++ll ++ + +fmqagf+l+e+g+vr+kn +n+++knl d++vv+++y++fGy l+ g+s         l++++ ++++e l+++   lf ++f
   WP_004745589.1  10 LWLLGCTCIALFMQAGFTLIETGSVRAKNSVNVAMKNLADFIVVTVAYIVFGYHLSQGDS--------LLSFETLAMDNEHLPKL---LFNVMF 92 
                      59999*****************************************************98........45567788999999888...6***** PP

                      HHHHHHHHHHHCTTTB-HHHHHHHHHHHHHHTHHHHHHHHHS.SSHHHHTT--.-SSSTTTTHHHHHHHHHHHHHHHS-STTTTG.TTTSS--G CS
  Ammonium_transp  95 aataitivsgavaerikfsayllfsallgtlvyppvahwvwgeggwlaklgvliDfAgstvVHlvggvagLaaalvlgkregrfe.gkeeaikg 187
                       +ta+tivsg vaer+  + y+  s++++ ++yp+ ++w+w+   wl++ g++ DfAg+t+VH+vgg++gL++++++g+r+grf+ +  ++i+ 
   WP_004745589.1  93 VTTAATIVSGCVAERMSYKGYIYTSFFIAVITYPIASYWTWNPYSWLNSSGFY-DFAGGTTVHVVGGMIGLVGTMIVGPRKGRFDsKSVREIPS 185
                      *****************************************************.*******************************99999**** PP

                      SBHHHHHHHHHHHHHHHHHHHHGGGSSSSHHHHH.HHHHHHHHHHHHHHHHHHHHHHHCSS-BHHHHHHHHHHHHHHHTTTTTTS-HHHHHHHH CS
  Ammonium_transp 188 hnlpfavlGtllLwfgWfgFNaGsaltankrararaavtTllAaaagaltallisrlkegkinvlglanGilAGlVAiTaacavvepwgAliiG 281
                      ++ ++++lG++l+ f W+gFN+Gs +t + r  + + ++Tl+ +a+++ t+l+   +   ++ v+ ++n++l+GlV +Ta+++++ ++  l++G
   WP_004745589.1 186 YSHTLVTLGVFLMLFAWLGFNGGSLYTFDLRVPK-ILFNTLVCGAIAGCTTLFWLHHY-RHVPVFVVLNSVLGGLVIVTAGADIMAMIDLLLLG 277
                      **********************************.***************99977766.89********************************* PP

                      HHHHHHHHHHHHHHHHHHTC-CTTHHHHHHHHHHHHHHHHHHHHTSHHHCSSHSSTTGGGT-C..HHHHHHHHHHHHHHHHHHHHHHHHHHHHH CS
  Ammonium_transp 282 lvAgvlsvlgvkklkeklkidDsldvvavHgvgGiwGllavgifaaekvvaskisggllsgeg..kqlvvqligilvilayafvvtlilllllk 373
                      + A++ ++lg   +  k+kidD++++++vH+++G++G+l  g+ +           g+++++   +ql +q+ g ++++++a++ ++++++ll 
   WP_004745589.1 278 MFASACVILGDR-MLIKAKIDDPVGAIPVHLFCGVVGTLYAGFKL-----------GWIENQDvvNQLLMQVTGLVMVAGWAALNAVTIFMLLR 359
                      **********96.9************************9988533...........47777777899*************************** PP

                      HHTTSB--HHHHHHSHHHHHHS CS
  Ammonium_transp 374 lllgLrvseeeeevglDvaehg 395
                      ++   rv+e eeevgl+v+ehg
   WP_004745589.1 360 KMHLDRVTEREEEVGLNVSEHG 381
                      ***99****************8 PP
'''

