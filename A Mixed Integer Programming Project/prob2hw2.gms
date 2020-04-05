set j       set of years                          /year0*year10/;
set jp(j)   set of years whos ends start a period /year0, year2, year5/;
set p       set of processes                      /p1*p3/;
set chem    set of chemicals                      /chem1*chem3/

parameter yr(j) year number
/ year0 0
  year1 1
  year2 2
  year3 3
  year4 4
  year5 5
  year6 6
  year7 7
  year8 8
  year9 9
  year10 10/;

parameter len(jp) the length of a period
/ year0 2
  year2 3
  year5 5 /;

parameter life(jp) the rest of the project life
/ year0 10
  year2 8
  year5 5 /;

Table alpha(jp,p) variable investment cost
         p1      p2      p3
  year0  1.38    2.72    1.76
  year2  1.56    3.22    2.34
  year5  1.78    4.60    2.84
;
Table beta(jp,p) fixed investment cost
         p1      p2      p3
  year0  85      73      110
  year2  95      82      125
  year5  112     102     148
;
Table opexf(j,p) operating expense factor
         p1      p2      p3
  year0  0.0     0.0     0.0
  year1  0.4     0.6     0.5
  year2  0.4     0.6     0.5
  year3  0.5     0.7     0.6
  year4  0.5     0.7     0.6
  year5  0.5     0.7     0.6
  year6  0.6     0.8     0.7
  year7  0.6     0.8     0.7
  year8  0.6     0.8     0.7
  year9  0.6     0.8     0.7
  year10 0.6     0.8     0.7
;
Table prices(j, chem) chemical prices
         chem1   chem2   chem3
  year0  0       0       0
  year1  4       9.6     26.20
  year2  4       9.6     26.20
  year3  5.24    11.52   29.20
  year4  5.24    11.52   29.20
  year5  5.24    11.52   29.20
  year6  7.32    13.52   35.20
  year7  7.32    13.52   35.20
  year8  7.32    13.52   35.20
  year9  7.32    13.52   35.20
  year10 7.32    13.52   35.20
;

parameter i      intrest rate per year  /0.10/;
parameter fsal   salvage value factor   /0.10/;
parameter fwc    working capital factor /0.15/;
parameter tax    tax rate               /0.45/;
parameter fpv(j) present value factor;

* calculate present value factor
loop(  j,  fpv(j) = 1.0/((1.0+i)**yr(j))   );


binary variable   y(jp,p)      build process p in period jp or not;
positive variable addCap(jp,p) additional capacity for process p in period jp;
positive variable cap(j, p)    total capacity of p in year j;
positive variable prflow(j,p)  product flows from process p in year j;
positive variable inflow(j,p)  reactant flow into p in year j;
positive variable pur1(j)      amount of chem 1 purchased year j;
positive variable pur2(j)      amount of chem 2 purchased year j;
positive variable sel3(j)      amount of chem 3 sold year j;
positive variable inv(jp)      amount invested in jp;
positive variable sell(j)      amount earned from the sale of chem3 in year j;
positive variable opex(j)      operating expenses for year j;
positive variable wc(jp)       working capital put in at begining of period jp;
positive variable buy1(j)      amount spent to purchase chem1 in year j;
positive variable buy2(j)      amount spent to purchace chem2 in year j;
positive variable income(j)    taxable imcome year j;
positive variable dep(j)       depriciation in year j;
variable npv                   net present value;

* set upper bound on purchases and sales
pur1.up(j) = 6$(yr(j) > yr('year0'))
             + 1.5$(yr(j) > yr('year2'))
             + 1.1$(yr(j) > yr('year5'));
pur2.up(j) = 20$(yr(j) > yr('year0'))
             + 5.5$(yr(j) > yr('year2'))
             + 4.5$(yr(j) > yr('year5'));
sel3.up(j) = 65$(yr(j) > yr('year0'))
             + 10$(yr(j) > yr('year2'))
             + 15$(yr(j) > yr('year5'));

* set upper bounds on investment amounts
inv.up('year0') = 200;
inv.up('year2') = 300;
inv.up('year5') = 400;

*
equations
   mb_p1(j)        mass balance process 1
   mb_p2(j)        mass balance process 2
   mb_p3(j)        mass balance process 3
   mb_p2p3(j)      mass balabce split between p2 and p3
   addcaplim(jp,p) limit on additional capacity added in a period
   p1cap(j)        eq for capacity of process 1 in year j
   p2cap(j)        eq for capacity of process 2 in year j
   p3cap(j)        eq for capacity of process 3 in year j



   inveq(jp)       eq for amount of investment in period jp
   npveq           eq for net present value
   opexeq(j)       eq for operating expenses
   wceq(jp)        eq for the amount of working capital put for period jp
   sel3eq(j)       eq for the amount of chem 3 sold in year j
   buy1eq(j)       eq for purchase costs of chem1 in year j
   buy2eq(j)       eq for purchase costs of chem2 in year j
   sell3eq(j)      eq for amount made from selling chem3
   prlim(j,p)      eq for capacity limit on product flows
   incomeeq(j)     eq for income in year j
   depeq(j)        eq for depriciation in year j
;


* Mass balances
mb_p1(j)  .. 1.11*prflow(j,'p1') =e= pur1(j);
mb_p2(j)  .. 1.22*prflow(j,'p2') =e= inflow(j,'p2');
mb_p3(j)  .. 1.05*prflow(j,'p3') =e= inflow(j,'p3');
mb_p2p3(j).. inflow(j,'p2') + inflow(j,'p3') =e= prflow(j,'p1') + pur2(j);

* Product capacity limits
prlim(j,p) .. prflow(j,p) =l= cap(j,p);
p1cap(j) .. cap(j,'p1') =e= sum(jp$(yr(j) > yr(jp)), addCap(jp,'p1'));
p2cap(j) .. cap(j,'p2') =e= sum(jp$(yr(j) > yr(jp)), addCap(jp,'p2'))
                            + 50$(yr(j) > 0);
p3cap(j) .. cap(j,'p3') =e= sum(jp$(yr(j) > yr(jp)), addCap(jp,'p3'));

* Investments and working capital
addcaplim(jp,p) .. addCap(jp,p) =l= 100*y(jp,p);
inveq(jp).. inv(jp) =e= sum(p, beta(jp,p)*y(jp,p))
                        + sum(p, alpha(jp,p)*addCap(jp,p));
wceq(jp) .. wc(jp)  =e= fwc*inv(jp);

* Buy and sell chemicals
sel3eq(j) .. sel3(j) =e= prflow(j,'p2')+prflow(j,'p3');
buy1eq(j) .. buy1(j) =e= pur1(j)*prices(j,'chem1');
buy2eq(j) .. buy2(j) =e= pur2(j)*prices(j,'chem2');
sell3eq(j).. sell(j) =e= sel3(j)*prices(j,'chem3');

* Income operating expences and depriciation
opexeq(j) .. opex(j) =e= sum(p, prflow(j,p)*opexf(j,p));
incomeeq(j) .. income(j) =e= sell(j) - opex(j)- buy1(j) - buy2(j);
depeq(j) .. dep(j) =e= sum(jp$(yr(j) > yr(jp)), (1.0-fsal)*inv(jp)/life(jp) );

* Calculate NPV
npveq .. npv =e= sum(j, dep(j)*fpv(j))*tax
                 + sum(j, income(j)*fpv(j))*(1-tax)
                 - sum(jp, fpv(jp)*(inv(jp) + wc(jp)))
                 + fpv('year10')*sum(jp, wc(jp) + fsal*inv(jp));

* Make sure I get the best answer
option optcr = 0;

* solve all equations maximize NPV
model mod /all/;
solve mod using mip maximizing npv;
