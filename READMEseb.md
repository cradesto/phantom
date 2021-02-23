соединение с leo

.ssh/config:

Host leo
   User seb
   Hostname 192.168.1.10
   Port 22
   ProxyJump gray
   HostKeyAlias leo
   CheckHostIP no

.alias:

a leo 'ssh -x seb@leo'
a leoy 'ssh -Y seb@leo'
a leox 'ssh -X seb@leo'
a toleo 'scp -pr \!* seb@leo:'
a frleo 'scp -pr seb@leo:\!* .'
a leoscr 'sshfs leo:/home/seb/ leo'

установка переменных окружения:

a setphantom 'setenv OMP_SCHEDULE "dynamic"; setenv OMP_STACKSIZE 512M; setenv SYSTEM gfortran'

a movie 'ffmpeg -i \!*_%04d.png -r 10 -vb 50M -bt 100M -pix_fmt yuv420p  movie.mp4'


look here
https://phantomsph.readthedocs.io/en/latest/running-first-calculation.html

/home/seb/prg/gitWork/phantom/scripts/writemake.sh  star > Makefile

mkphst  star
mkpoly polytrope

вместо star может быть blast, polytrope, nsmerger, etc.

before make setup do
 setphantom

before running phantom to relax  star insert
idam = 1 (or 2?)
damp =       0.030    ! artificial damping of velocities (if on, v=0 initially)


Если выбрать единицу массы Msun, а длины 1 км, задать 0.3 + 1.4 Msun, то дальше по умолчанию:

---------- binary parameters -----------
  primary mass     :   0.300
  secondary mass   :    1.40
  mass ratio m2/m1 :    4.67
  reduced mass     :   0.247
  semi-major axis  :    100.
  period           :   0.482E+04
  eccentricity     :    0.00
  pericentre       :    100.
  apocentre        :    100.
  angular momentum :   3.221
  mean ang. speed  :  0.1304E-02
  Omega_0 (prim)   :  0.1304E-02
  Omega_0 (second) :  0.1304E-02
  R_accretion (1)  :   1.000
  R_accretion (2)  :   1.000
  Roche lobe  (1)  :   25.66
  Roche lobe  (2)  :   51.48
  ----------------------------------------

Беру gamma=5/3
alpha_min = 0.1
Rkill = 60 km
R2=15 (M2=1.4)



Andrey Yudin

Я бы запустил со следующими параметрами (хотя УрС и не тот):
m1=1.4, R1=11 km
m2=0.4, R2=13 km
a=46.4 km - расстояние между компаньонами. Тогда m2 чуть переполняет свою
полость Роша.

Беру
m1=1.2, R1=11 km
m2=0.4, R2=13 km -- чтобы делить на 3 по числу точек
a=48.3497 km - расстояние между компаньонами. Тогда m2 чуть переполняет свою
полость Роша.
gamma=5/3, a~50 быстрый разрыв
а~60 разрыв помедленнее


1) gamma=2
2) gamma=5/3

для расчёта релаксации нужно alpha=0.1, потом в binary можно 0

binary has wrong description at make moddump -- оно даёт sink вместо реальной 2й звезды, но это часто удобно.
-- такой режим по умолчанию для mkphst (star)
make moddump
is in fact
make moddump MODFILE=moddump_binary.f90
which does note work for mkpoly regime

it must be

 make moddump MODFILE=moddump_binarystar.f90
-- такой режим по умолчанию для mkpoly (polytrope)

for magnetic field do
 make moddump MODFILE=moddump_binarystar.f90 MHD=yes

then

 phantommoddump star_00100 binary_00000.tmp 0.0
или
 phantommoddump nstar_00100 binary_00000.tmp 0.0
etc

Code units of mag flux density =    1.8213E+020

Initial field strength (Bzero) =  -1.0000E+13 G ( -5.4905E-08 in code units)

 --- code units ---

     Mass:  1.989E+33 g       Length:  1.000E+05 cm    Time:  2.745E-06 s
  Density:  1.989E+18 g/cm^3  Energy:  2.640E+54 erg   En/m:  1.327E+21 erg/g
 Velocity:  3.643E+10 cm/s    Bfield:  1.821E+20 G
        G:  1.000E+00              c:  8.229E-01       mu_0:  1.000E+00

 Column Density:  1.989E+23 g/cm^2

масса измеряется в Msun, а длина в км.
Time unit = sqrt(Length^3/(G_N*Msun))


For TDE (tde2) with M_BH=1e6 Msun we have M_u=Msun, L_u=Rsun:
Mstar=M_*=25 Msun  e=0 - эксцентиситет
orbit 300 Rsun Raccr for BH 13 Rsun

--- code units ---

     Mass:  1.989E+33 g       Length:  6.960E+10 cm    Time:  1.594E+03 s
  Density:  5.901E+00 g/cm^3  Energy:  3.793E+48 erg   En/m:  1.907E+15 erg/g
 Velocity:  4.367E+07 cm/s    Bfield:  3.760E+08 G
        G:  1.000E+00              c:  6.865E+02       mu_0:  1.000E+00

Для M_BH=1e6 Msun и d=300Rsun=2.1e13 cm= 1.4 AU по 3му закону Кеплера:

T^2=d^3/(M_{BH}+M_*), если T в годах, а d в AU, T=1.6e-3 yr

tde3 M_BH=1e6 Msun M_*=25 Msun e=0.15
orbit 300 Rsun Raccr for BH 13 Rsun

Enter accretion radius for the companion in code units ([0.000:], default=0.000): 13
Do you want a corotating frame with a corotating binary? (default=no):
 reset CofM: (-1.49E-13 -5.59E-14 -4.38E-14 ) -> (-3.52E-17  6.18E-17 -1.70E-17 )
  Got       100000      100000  after deleting accreted particles

  ---------- binary parameters -----------
  primary mass     :    25.0
  secondary mass   :   0.100E+07
  mass ratio m2/m1 :   0.400E+05
  reduced mass     :    25.0
  semi-major axis  :    300.
  period           :    32.6
  eccentricity     :   0.150
  pericentre       :    255.
  apocentre        :    345.
  angular momentum :  0.4281E+06
  mean ang. speed  :  0.1655
  Omega_0 (prim)   :  0.1439
  Omega_0 (second) :  0.1439
  R_accretion (1)  :   0.000
  R_accretion (2)  :   13.00
  Roche lobe  (1)  :   4.285
  Roche lobe  (2)  :   243.8
  ----------------------------------------

tde4 all the same but e=0.5:
 ---------- binary parameters -----------
  primary mass     :    25.0
  secondary mass   :   0.100E+07
  mass ratio m2/m1 :   0.400E+05
  reduced mass     :    25.0
  semi-major axis  :    300.
  period           :    32.6
  eccentricity     :   0.500
  pericentre       :    150.
  apocentre        :    450.
  angular momentum :  0.3750E+06
  mean ang. speed  :  0.1111
  Omega_0 (prim)   :  0.7407E-01
  Omega_0 (second) :  0.7407E-01
  R_accretion (1)  :   0.000
  R_accretion (2)  :   13.00
  Roche lobe  (1)  :   4.285
  Roche lobe  (2)  :   243.8
  ----------------------------------------

 tde5 e=0.7

 ---------- binary parameters -----------
  primary mass     :    25.0
  secondary mass   :   0.100E+07
  mass ratio m2/m1 :   0.400E+05
  reduced mass     :    25.0
  semi-major axis  :    300.
  period           :    32.6
  eccentricity     :   0.700
  pericentre       :    90.0
  apocentre        :    510.
  angular momentum :  0.3092E+06
  mean ang. speed  :  0.8085E-01
  Omega_0 (prim)   :  0.4756E-01
  Omega_0 (second) :  0.4756E-01
  R_accretion (1)  :   0.000
  R_accretion (2)  :   13.00
  Roche lobe  (1)  :   4.285
  Roche lobe  (2)  :   243.8
  ---------------------------------------

tde6 e=0.9

 ---------- binary parameters -----------
  primary mass     :    25.0
  secondary mass   :   0.100E+07
  mass ratio m2/m1 :   0.400E+05
  reduced mass     :    25.0
  semi-major axis  :    300.
  period           :    32.6
  eccentricity     :   0.900
  pericentre       :    30.0
  apocentre        :    570.
  angular momentum :  0.1887E+06
  mean ang. speed  :  0.4415E-01
  Omega_0 (prim)   :  0.2324E-01
  Omega_0 (second) :  0.2324E-01
  R_accretion (1)  :   0.000
  R_accretion (2)  :   13.00
  Roche lobe  (1)  :   4.285
  Roche lobe  (2)  :   243.8
  ----------------------------------------


 Полный момент импульса сохраняется?

 В описании они пишут:
 Typically with individual timesteps, one should expect
 energy conservation to  E/E ∼ 10^{−3} and linear and angular
 momentum conservation to ∼10^{−6}} with default code settings.
 The code execution is aborted if conservation errors exceed
 10%.

Я нашёл в файле .ev -- монитор эволюции. Я в нём давно искал, но
колонка 170 была вне поля зрения, называется angmom.

в начале было
2.8386115

в конце
2.8384699

-- неплохо сохраняется.



 plot
 ssplash binary_00000.tmp

then run
 это не нужно: ./phantomsetup binary.in
 ./phantom binary.in

unit km actually is 10 km

Пример
 mkdir binNS
 cd binNS/
 locate writemake.sh
 /home/seb/prg/gitWork/phantom/scripts/writemake.sh star > Makefile
делаем первую нейтронную звезду
 make setup
 j star.setup
 ./phantomsetup star
 j star.in
 make
 phantom star.in
 -- пускаем phantom, чтобы срелаксировать её до равновесия, можно смотреть:
      ssplash star_000*
      make moddump MODFILE=moddump_binarystar.f90

переименуем выдачи star_* в nsA_#1, а star.setup и star.in сотрём.
Снова пустим ./phantomsetup star
j star.in
phantom star.in

phantommoddump star_00030 binNS_00000.tmp 0.0

a setphantom 'setenv OMP_SCHEDULE "dynamic";setenv OMP_STACKSIZE 512M;setenv SYSTEM gfortran'
a mkphst '/home/seb/prg/gitWork/phantom/scripts/writemake.sh star > Makefile'
a mkpoly '/home/seb/prg/gitWork/phantom/scripts/writemake.sh polytrope > Makefile'


Для нейтр.звёзд урсос Павла Хензела:
 cd composit/
 make clean
 del *
 mkpoly
 make setup
 phantomsetup nsB
 6) Piecewise polytrope
-- выбираем broken polytrop, Eos 9, а там SLy:
 2=SLy

 nsA:
 Radius              =  1.25760E+06 cm     =  1.25760E+01 km
 Mass                =  2.78474E+33 g      =  1.40000E+00 M_sun
 N                   =       100000

 phantomsetup nsB
 -- выбираем relax
 Radius              =  9.79200E+05 cm     =  9.79200E+00 km
 Mass                =  7.95640E+32 g      =  4.00000E-01 M_sun
 N                   =       100000

 make
 phantom nsB.in


 ns01:
 Radius              =  7.48800E+05 cm     =  7.48800E+00 km
 Mass                =  1.98910E+32 g      =  1.00000E-01 M_sun
 N                   =       100000 (10000 - то же)

1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)


 ns01A APR3
 Radius              =  7.68000E+05 cm     =  7.68000E+00 km
 Mass                =  1.98910E+32 g      =  1.00000E-01 M_sun

 ns01C MS1
 Radius              =  8.51400E+05 cm     =  8.51400E+00 km
 Mass                =  1.98910E+32 g      =  1.00000E-01 M_sun

 ns01C ENG
 Radius              =  6.91200E+05 cm     =  6.91200E+00 km
 Mass                =  1.98910E+32 g      =  1.00000E-01 M_sun

For n=2 and requested R=1e4 final
 Radius              =  1.11941E+07 cm     =  1.11941E+02 km
 Mass                =  7.95640E+32 g      =  4.00000E-01 M_sun
 rho_central         =  4.45240E+11 g/cm^3
 N                   =       200000

Для меньшей массы:
gamma               =      2.00000
 polyk               =   7980.35108
setting isothermal sound speed^2 (polyk) =    7980.3510812959958       gamma =    2.00

 Radius              =  1.11941E+07 cm     =  1.11941E+02 km
 Mass                =  3.97820E+32 g      =  2.00000E-01 M_sun

phantommoddump ns02n1_00100 d400m04m02_00000 0.0
using ns04n1_00100 as a second star

 q=   0.50000     RL1=   0.5707515723
 Ry=   0.4399444148      Rz=   0.4142650888
Тогда расстояние должно быть 270.2 км при заполнении полости Роша малой массой.
1st run with 300 km
2nd run with 320 km
3rd run with 400 km
4th run with 360 km

магнитим меньшую звезду:
phantommoddump ns02n1_00100 magn02n1_00000.tmp 0.0

phantommoddump magn02n1_00000  magnd360m04m02_00000.tmp 0.0

-- эта команда сажает магн поле из звезды magn02n1 на
звезду 0.4 Msun.

приходится делать так:

phantommoddump ns04n1_00100 magnd360m04m02_00000.tmp 0.0

phantommoddump ns04n1_00100 mag1e4d360m04m02_00000.tmp 0.0


1) B=1e13 Gs no effect

2) B=1e15
