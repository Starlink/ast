/*
*  Purpose:
*     Check the Moc class against the CDS HEALPix library.
*
*     The expected values in the tables below were computed with the Rust
*     libraries published by the CDS: cdshealpix 0.10.1 as used by cds-moc
*     0.19.3 (https://github.com/cds-astro/cds-moc-rust), the library that
*     MOCPy is built on.  They check that AST and the CDS libraries agree on:
*
*     - the HEALPix cell that holds a sky position, at orders from 0 to 27.
*       Many of the positions are on the boundary of a cell, and on a
*       meridian or at a pole where the projections that the Moc class uses
*       are cut, since those are the positions that are easy to get wrong;
*     - the cells of a Moc made from a Circle, including Circles centred on
*       those positions, Circles that cover more than half the sky, and Mocs
*       with a maximum order of zero;
*     - the cells of a Moc made from a Polygon;
*     - the string and JSON serialisations of a MOC.
*
*  Notes:
*     - Public C API only, with astWatch/astOK for error checking.
*     - A position on the boundary of a cell is in more than one cell, and a
*       library may put it in any of them.  Each entry in the table of
*       positions lists the cells that cdshealpix gives for the position and
*       for positions displaced from it by a small fraction of a cell in
*       eight directions.  AST must put the position in one of them: a Moc
*       made of those cells must contain it and a Moc made of all the other
*       cells must not.
*     - A Moc made from a Region holds the cells at the maximum order whose
*       centre is inside the Region, so long as MinOrder is equal to MaxOrder
*       (otherwise parts of a Region smaller than a cell at MinOrder can be
*       missed, as astAddRegion documents).  The expected cells of each
*       Circle and Polygon are those whose centre, as given by cdshealpix, is
*       inside it: within the radius of the Circle by angular distance, or on
*       the inner side of the great circle through each edge of the Polygon.
*       The cdshealpix functions that return the cells covering a cone or
*       polygon were not used for this, since they return cells that overlap
*       it, or are approximate.
*     - Cases in which the centre of any cell is within 1E-6 degrees of the
*       edge of the Region were left out, so that no expected value depends
*       on rounding.
*     - The vertices of each Polygon are in the order that makes AST take the
*       smaller of the two areas they bound as the inside.
*/
#include "ast.h"
#include "ast_err.h"

#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/* A sky position, in degrees, and the cells at the given order that hold
   it. */
typedef struct {
   int order;
   double lon;
   double lat;
   int ncell;
   int64_t cell[ 4 ];
} Position;

/* A cone, in degrees, and the string MOC of the cells at the given order
   whose centre is inside it. */
typedef struct {
   int order;
   double lon;
   double lat;
   double radius;
   const char *moc;
} Cone;

/* A polygon, with vertices given in degrees as longitude then latitude, and
   the string MOC of the cells at the given order whose centre is inside
   it. */
typedef struct {
   int order;
   int nvert;
   double vert[ 10 ];
   const char *moc;
} Poly;

/* The reference values. */
static const Position positions[] = {
   { 0, 0.0, 0.0, 1, { INT64_C(4) } },
   { 0, 90.0, 60.0, 2, { INT64_C(0), INT64_C(1) } },
   { 0, 180.0, -75.0, 2, { INT64_C(9), INT64_C(10) } },
   { 0, 270.0, 41.810314895778596, 3,
     { INT64_C(2), INT64_C(3), INT64_C(7) } },
   { 0, 45.0, 90.0, 1, { INT64_C(0) } },
   { 0, 225.0, -90.0, 1, { INT64_C(10) } },
   { 0, 0.0, 41.810314895778596, 3, { INT64_C(0), INT64_C(3), INT64_C(4) } },
   { 0, 90.0, 89.9999999, 2, { INT64_C(0), INT64_C(1) } },
   { 0, 359.9999999, 5.0, 1, { INT64_C(4) } },
   { 0, 135.0, -60.0, 1, { INT64_C(9) } },
   { 0, 0.0, -30.0, 1, { INT64_C(4) } },
   { 0, 180.0, 0.0, 1, { INT64_C(6) } },
   { 0, 270.0, -41.810314895778596, 3,
     { INT64_C(7), INT64_C(10), INT64_C(11) } },
   { 0, 315.0, 90.0, 1, { INT64_C(3) } },
   { 0, 90.0, 0.0, 1, { INT64_C(5) } },
   { 0, 0.0, 60.0, 2, { INT64_C(0), INT64_C(3) } },
   { 0, 22.5, 19.47122063449069, 2, { INT64_C(0), INT64_C(4) } },
   { 0, 123.456, -54.321, 1, { INT64_C(9) } },
   { 0, 301.2, 77.7, 1, { INT64_C(3) } },
   { 0, 67.5, -19.47122063449069, 2, { INT64_C(5), INT64_C(8) } },
   { 1, 0.0, 0.0, 4, { INT64_C(16), INT64_C(17), INT64_C(18), INT64_C(19) } },
   { 1, 90.0, 60.0, 2, { INT64_C(1), INT64_C(6) } },
   { 1, 180.0, -75.0, 2, { INT64_C(36), INT64_C(40) } },
   { 1, 270.0, 41.810314895778596, 3,
     { INT64_C(9), INT64_C(14), INT64_C(31) } },
   { 1, 45.0, 90.0, 1, { INT64_C(3) } },
   { 1, 225.0, -90.0, 1, { INT64_C(40) } },
   { 1, 0.0, 41.810314895778596, 3,
     { INT64_C(2), INT64_C(13), INT64_C(19) } },
   { 1, 90.0, 89.9999999, 2, { INT64_C(3), INT64_C(7) } },
   { 1, 359.9999999, 5.0, 1, { INT64_C(19) } },
   { 1, 135.0, -60.0, 1, { INT64_C(36) } },
   { 1, 0.0, -30.0, 1, { INT64_C(16) } },
   { 1, 180.0, 0.0, 4,
     { INT64_C(24), INT64_C(25), INT64_C(26), INT64_C(27) } },
   { 1, 270.0, -41.810314895778596, 3,
     { INT64_C(28), INT64_C(41), INT64_C(46) } },
   { 1, 315.0, 90.0, 1, { INT64_C(15) } },
   { 1, 90.0, 0.0, 4,
     { INT64_C(20), INT64_C(21), INT64_C(22), INT64_C(23) } },
   { 1, 0.0, 60.0, 2, { INT64_C(2), INT64_C(13) } },
   { 1, 22.5, 19.47122063449069, 4,
     { INT64_C(0), INT64_C(2), INT64_C(17), INT64_C(19) } },
   { 1, 123.456, -54.321, 1, { INT64_C(36) } },
   { 1, 301.2, 77.7, 1, { INT64_C(15) } },
   { 1, 67.5, -19.47122063449069, 4,
     { INT64_C(20), INT64_C(22), INT64_C(33), INT64_C(35) } },
   { 3, 0.0, 0.0, 4,
     { INT64_C(271), INT64_C(282), INT64_C(293), INT64_C(304) } },
   { 3, 90.0, 60.0, 2, { INT64_C(29), INT64_C(110) } },
   { 3, 180.0, -75.0, 2, { INT64_C(580), INT64_C(648) } },
   { 3, 270.0, 41.810314895778596, 3,
     { INT64_C(149), INT64_C(234), INT64_C(511) } },
   { 3, 45.0, 90.0, 1, { INT64_C(63) } },
   { 3, 225.0, -90.0, 1, { INT64_C(640) } },
   { 3, 0.0, 41.810314895778596, 3,
     { INT64_C(42), INT64_C(213), INT64_C(319) } },
   { 3, 90.0, 89.9999999, 2, { INT64_C(63), INT64_C(127) } },
   { 3, 359.9999999, 5.0, 1, { INT64_C(304) } },
   { 3, 135.0, -60.0, 1, { INT64_C(588) } },
   { 3, 0.0, -30.0, 4,
     { INT64_C(256), INT64_C(257), INT64_C(258), INT64_C(259) } },
   { 3, 180.0, 0.0, 4,
     { INT64_C(399), INT64_C(410), INT64_C(421), INT64_C(432) } },
   { 3, 270.0, -41.810314895778596, 3,
     { INT64_C(448), INT64_C(661), INT64_C(746) } },
   { 3, 315.0, 90.0, 1, { INT64_C(255) } },
   { 3, 90.0, 0.0, 4,
     { INT64_C(335), INT64_C(346), INT64_C(357), INT64_C(368) } },
   { 3, 0.0, 60.0, 2, { INT64_C(46), INT64_C(221) } },
   { 3, 22.5, 19.47122063449069, 4,
     { INT64_C(10), INT64_C(32), INT64_C(287), INT64_C(309) } },
   { 3, 123.456, -54.321, 1, { INT64_C(590) } },
   { 3, 301.2, 77.7, 1, { INT64_C(254) } },
   { 3, 67.5, -19.47122063449069, 4,
     { INT64_C(330), INT64_C(352), INT64_C(543), INT64_C(565) } },
   { 8, 0.0, 0.0, 4,
     { INT64_C(278527), INT64_C(289450), INT64_C(300373), INT64_C(311296) } },
   { 8, 90.0, 60.0, 2, { INT64_C(30711), INT64_C(113659) } },
   { 8, 180.0, -75.0, 2, { INT64_C(594177), INT64_C(664066) } },
   { 8, 270.0, 41.810314895778596, 3,
     { INT64_C(152917), INT64_C(240298), INT64_C(524287) } },
   { 8, 45.0, 90.0, 1, { INT64_C(65535) } },
   { 8, 225.0, -90.0, 1, { INT64_C(655360) } },
   { 8, 0.0, 41.810314895778596, 3,
     { INT64_C(43690), INT64_C(218453), INT64_C(327679) } },
   { 8, 90.0, 89.9999999, 2, { INT64_C(65535), INT64_C(131071) } },
   { 8, 359.9999999, 5.0, 1, { INT64_C(312064) } },
   { 8, 135.0, -60.0, 1, { INT64_C(602883) } },
   { 8, 0.0, -30.0, 4,
     { INT64_C(263167), INT64_C(263850), INT64_C(264533), INT64_C(265216) } },
   { 8, 180.0, 0.0, 4,
     { INT64_C(409599), INT64_C(420522), INT64_C(431445), INT64_C(442368) } },
   { 8, 270.0, -41.810314895778596, 3,
     { INT64_C(458752), INT64_C(677205), INT64_C(764586) } },
   { 8, 315.0, 90.0, 1, { INT64_C(262143) } },
   { 8, 90.0, 0.0, 4,
     { INT64_C(344063), INT64_C(354986), INT64_C(365909), INT64_C(376832) } },
   { 8, 0.0, 60.0, 2, { INT64_C(48123), INT64_C(227319) } },
   { 8, 22.5, 19.47122063449069, 4,
     { INT64_C(10922), INT64_C(32768), INT64_C(294911), INT64_C(316757) } },
   { 8, 123.456, -54.321, 1, { INT64_C(604821) } },
   { 8, 301.2, 77.7, 1, { INT64_C(260496) } },
   { 8, 67.5, -19.47122063449069, 4,
     { INT64_C(338602), INT64_C(360448), INT64_C(557055), INT64_C(578901) } },
   { 14, 0.0, 0.0, 4,
     { INT64_C(1140850687), INT64_C(1185589930), INT64_C(1230329173),
       INT64_C(1275068416) } },
   { 14, 90.0, 60.0, 2, { INT64_C(125795829), INT64_C(465551098) } },
   { 14, 180.0, -75.0, 2, { INT64_C(2433750292), INT64_C(2720016936) } },
   { 14, 270.0, 41.810314895778596, 3,
     { INT64_C(626349397), INT64_C(984263338), INT64_C(2147483647) } },
   { 14, 45.0, 90.0, 1, { INT64_C(268435455) } },
   { 14, 225.0, -90.0, 1, { INT64_C(2684354560) } },
   { 14, 0.0, 41.810314895778596, 3,
     { INT64_C(178956970), INT64_C(894784853), INT64_C(1342177279) } },
   { 14, 90.0, 89.9999999, 2, { INT64_C(268435455), INT64_C(536870911) } },
   { 14, 359.9999999, 5.0, 1, { INT64_C(1278217468) } },
   { 14, 135.0, -60.0, 1, { INT64_C(2469408963) } },
   { 14, 0.0, -30.0, 4,
     { INT64_C(1077936127), INT64_C(1080732330), INT64_C(1083528533),
       INT64_C(1086324736) } },
   { 14, 180.0, 0.0, 4,
     { INT64_C(1677721599), INT64_C(1722460842), INT64_C(1767200085),
       INT64_C(1811939328) } },
   { 14, 270.0, -41.810314895778596, 3,
     { INT64_C(1879048192), INT64_C(2773833045), INT64_C(3131746986) } },
   { 14, 315.0, 90.0, 1, { INT64_C(1073741823) } },
   { 14, 90.0, 0.0, 4,
     { INT64_C(1409286143), INT64_C(1454025386), INT64_C(1498764629),
       INT64_C(1543503872) } },
   { 14, 0.0, 60.0, 2, { INT64_C(197115642), INT64_C(931102197) } },
   { 14, 22.5, 19.47122063449069, 4,
     { INT64_C(44739242), INT64_C(134217728), INT64_C(1207959551),
       INT64_C(1297438037) } },
   { 14, 123.456, -54.321, 1, { INT64_C(2477349348) } },
   { 14, 301.2, 77.7, 1, { INT64_C(1066993847) } },
   { 14, 67.5, -19.47122063449069, 4,
     { INT64_C(1386916522), INT64_C(1476395008), INT64_C(2281701375),
       INT64_C(2371179861) } },
   { 27, 0.0, 0.0, 4,
     { INT64_C(76561193665298431), INT64_C(79563593416878762),
       INT64_C(82565993168459093), INT64_C(85568392920039424) } },
   { 27, 90.0, 60.0, 2,
     { INT64_C(8442015247097821), INT64_C(31242605387771886) } },
   { 27, 180.0, -75.0, 2,
     { INT64_C(163326217360262416), INT64_C(182537246644668960) } },
   { 27, 270.0, 41.810314895778596, 3,
     { INT64_C(42033596522124629), INT64_C(66052794534767274),
       INT64_C(144115188075855871) } },
   { 27, 45.0, 90.0, 1, { INT64_C(18014398509481983) } },
   { 27, 225.0, -90.0, 1, { INT64_C(180143985094819840) } },
   { 27, 0.0, 41.810314895778596, 3,
     { INT64_C(12009599006321322), INT64_C(60047995031606613),
       INT64_C(90071992547409919) } },
   { 27, 90.0, 89.9999999, 2,
     { INT64_C(18014398509481983), INT64_C(36028797018963967) } },
   { 27, 359.9999999, 5.0, 1, { INT64_C(85779722289479872) } },
   { 27, 135.0, -60.0, 3,
     { INT64_C(165719230308731916), INT64_C(165719230308731917),
       INT64_C(165719230308731918) } },
   { 27, 0.0, -30.0, 4,
     { INT64_C(72339069014638591), INT64_C(72526718999112362),
       INT64_C(72714368983586133), INT64_C(72902018968059904) } },
   { 27, 180.0, 0.0, 4,
     { INT64_C(112589990684262399), INT64_C(115592390435842730),
       INT64_C(118594790187423061), INT64_C(121597189939003392) } },
   { 27, 270.0, -41.810314895778596, 3,
     { INT64_C(126100789566373888), INT64_C(186148784597980501),
       INT64_C(210167982610623146) } },
   { 27, 315.0, 90.0, 1, { INT64_C(72057594037927935) } },
   { 27, 90.0, 0.0, 4,
     { INT64_C(94575592174780415), INT64_C(97577991926360746),
       INT64_C(100580391677941077), INT64_C(103582791429521408) } },
   { 27, 0.0, 60.0, 2,
     { INT64_C(13228206878289902), INT64_C(62485210775543773) } },
   { 27, 22.5, 19.47122063449069, 4,
     { INT64_C(3002399751580330), INT64_C(9007199254740992),
       INT64_C(81064793292668927), INT64_C(87069592795829589) } },
   { 27, 123.456, -54.321, 1, { INT64_C(166252100488445844) } },
   { 27, 301.2, 77.7, 1, { INT64_C(71604745000766622) } },
   { 27, 67.5, -19.47122063449069, 4,
     { INT64_C(93074392298990250), INT64_C(99079191802150912),
       INT64_C(153122387330596863), INT64_C(159127186833757525) } },
};

static const Cone cones[] = {
   { 0, 0.0, 0.0, 40.0,
     "0/4" },
   { 0, 0.0, 0.0, 100.0,
     "0/0 3-5 7-8 11" },
   { 1, 0.0, 0.0, 40.0,
     "0/4 1/" },
   { 2, 0.0, 0.0, 40.0,
     "0/4 2/2 8 49 52 139 142 183 189" },
   { 3, 0.0, 0.0, 40.0,
     "0/4 3/8 10 32 34 40 196-197 208-209 212 555 558-559 570-571 727 733 "
     "735 757 759" },
   { 4, 0.0, 0.0, 12.0,
     "3/271 282 293 304 4/1078-1079 1081 1083 1122-1123 1132 1134 1169 "
     "1171 1180-1181 1220 1222 1224-1225" },
   { 5, 0.0, 0.0, 3.0,
     "5/4349-4351 4522 4693 4864-4866" },
   { 2, 0.0, 0.0, 100.0,
     "0/0 3-4 8 11 1/22 29 2/26-27 30-31 37 39 45 47 80 82-83 92 94-95 "
     "112-113 115 124-125 127 144 146 152 154 160-161 164-165" },
   { 2, 0.0, 0.0, 89.9,
     "0/0 3-4 8 11 1/22 29 2/82 94 113 125" },
   { 1, 0.0, 0.0, 120.0,
     "0/0 3-5 7-8 11 1/6-7 9 11 36 38 40-41" },
   { 0, 0.0, 90.0, 100.0,
     "0/0-7" },
   { 1, 0.0, 90.0, 40.0,
     "1/3 7 11 15" },
   { 2, 0.0, 90.0, 40.0,
     "1/3 7 11 15 2/7 11 23 27 39 43 55 59" },
   { 3, 0.0, 90.0, 40.0,
     "2/13-15 29-31 45-47 61-63 3/29-31 45-47 49-51 93-95 109-111 113-115 "
     "157-159 173-175 177-179 221-223 237-239 241-243" },
   { 4, 0.0, 90.0, 12.0,
     "3/63 127 191 255 4/245-247 249-251 501-503 505-507 757-759 761-763 "
     "1013-1015 1017-1019" },
   { 5, 0.0, 90.0, 3.0,
     "5/1021-1023 2045-2047 3069-3071 4093-4095" },
   { 2, 0.0, 90.0, 100.0,
     "0/0-3 1/17-19 21-23 25-27 29-31 2/67 83 99 115 143 159 175 191" },
   { 2, 0.0, 90.0, 150.0,
     "0/0-7 1/33-35 37-39 41-43 45-47 2/131 147 163 179" },
   { 2, 0.0, 90.0, 89.9,
     "0/0-3 1/19 23 27 31 2/71 75 87 91 103 107 119 123" },
   { 1, 0.0, 90.0, 120.0,
     "0/0-7 1/35 39 43 47" },
   { 0, 0.0, -90.0, 100.0,
     "0/4-11" },
   { 1, 0.0, -90.0, 40.0,
     "1/32 36 40 44" },
   { 2, 0.0, -90.0, 40.0,
     "1/32 36 40 44 2/132 136 148 152 164 168 180 184" },
   { 3, 0.0, -90.0, 40.0,
     "2/128-130 144-146 160-162 176-178 3/524-526 528-530 544-546 588-590 "
     "592-594 608-610 652-654 656-658 672-674 716-718 720-722 736-738" },
   { 4, 0.0, -90.0, 12.0,
     "3/512 576 640 704 4/2052-2054 2056-2058 2308-2310 2312-2314 "
     "2564-2566 2568-2570 2820-2822 2824-2826" },
   { 5, 0.0, -90.0, 3.0,
     "5/8192-8194 9216-9218 10240-10242 11264-11266" },
   { 2, 0.0, -90.0, 100.0,
     "0/8-11 1/16-18 20-22 24-26 28-30 2/0 16 32 48 76 92 108 124" },
   { 2, 0.0, -90.0, 150.0,
     "0/4-11 1/0-2 4-6 8-10 12-14 2/12 28 44 60" },
   { 2, 0.0, -90.0, 89.9,
     "0/8-11 1/16 20 24 28 2/68 72 84 88 100 104 116 120" },
   { 1, 0.0, -90.0, 120.0,
     "0/4-11 1/0 4 8 12" },
   { 0, 90.0, 89.0, 100.0,
     "0/0-7" },
   { 1, 90.0, 89.0, 40.0,
     "1/3 7 11 15" },
   { 2, 90.0, 89.0, 40.0,
     "1/3 7 11 15 2/7 11 23 27 39 43 55 59" },
   { 3, 90.0, 89.0, 40.0,
     "2/13-15 29-31 45-47 61-63 3/29-31 45-47 49-51 93-95 109-111 113-115 "
     "157-159 173-175 177-179 221-223 237-239 241-243" },
   { 4, 90.0, 89.0, 12.0,
     "3/63 127 191 255 4/245-247 249-251 501-503 505-507 759 762-763 1013 "
     "1015 1019" },
   { 5, 90.0, 89.0, 3.0,
     "5/1021-1023 2045-2047 3071 4095" },
   { 2, 90.0, 89.0, 100.0,
     "0/0-3 1/17-19 21-23 25-27 31 2/67 83 99 117-119 121-123 143 159" },
   { 2, 90.0, 89.0, 150.0,
     "0/0-7 1/33-35 37-39 41-43 45-47 2/131 147 163 179" },
   { 2, 90.0, 89.0, 89.9,
     "0/0-3 1/19 23 27 31 2/69-71 75 85-87 89-91 103 105-107 119 123" },
   { 1, 90.0, 89.0, 120.0,
     "0/0-7 1/35 39 43 47" },
   { 0, 270.0, -89.0, 100.0,
     "0/4-11" },
   { 1, 270.0, -89.0, 40.0,
     "1/32 36 40 44" },
   { 2, 270.0, -89.0, 40.0,
     "1/32 36 40 44 2/132 136 148 152 164 168 180 184" },
   { 3, 270.0, -89.0, 40.0,
     "2/128-130 144-146 160-162 176-178 3/524-526 528-530 544-546 588-590 "
     "592-594 608-610 652-654 656-658 672-674 716-718 720-722 736-738" },
   { 4, 270.0, -89.0, 12.0,
     "3/512 576 640 704 4/2052 2056 2058 2308-2309 2312 2564-2566 "
     "2568-2570 2820-2822 2824-2826" },
   { 5, 270.0, -89.0, 3.0,
     "5/8192 9216 10240-10242 11264-11266" },
   { 2, 270.0, -89.0, 100.0,
     "0/8-11 1/16-18 20 24-26 28-30 2/32 48 76 84-86 88-90 108 124" },
   { 2, 270.0, -89.0, 150.0,
     "0/4-11 1/0-2 4-6 8-10 12-14 2/12 28 44 60" },
   { 2, 270.0, -89.0, 89.9,
     "0/8-11 1/16 20 24 28 2/68 72-74 84 88 100-102 104 116-118 120-122" },
   { 1, 270.0, -89.0, 120.0,
     "0/4-11 1/0 4 8 12" },
   { 0, 180.0, 0.0, 40.0,
     "0/6" },
   { 0, 90.0, 60.0, 40.0,
     "0/0-1" },
   { 0, 45.0, 41.810314895778596, 40.0,
     "0/0" },
   { 2, 45.0, 41.810314895778596, 150.0,
     "0/0-9 11 2/160 175" },
   { 0, 0.0, 30.0, 40.0,
     "0/0 3-4" },
   { 1, 42.883, 0.288, 150.0,
     "0/0-1 3-5 8-9 11 1/9-11 24 26-29 31 40-42" },
   { 4, 279.566, -1.042, 170.0,
     "0/0-4 6-11 1/22 2/80-82 84-85 93-95 3/332 334 348-349 351 370-371 "
     "4/1332-1334 1340 1342-1343 1381 1401 1403 1474-1475 1478-1479" },
   { 1, 211.148, -31.341, 150.0,
     "0/1-11 1/1 3" },
   { 1, 289.849, 87.875, 150.0,
     "0/0-7 1/33-35 37-39 41-43 45-47" },
   { 4, 312.201, -8.576, 170.0,
     "0/0 2-4 6-11 1/5-7 20 22-23 2/17 19 84 86-87 3/74-75 340 342 4/261 "
     "290-291 293-295 1364 1366 1372 1374" },
   { 4, 143.333, 2.826, 170.0,
     "0/0-2 5-10 1/13-17 19 44-46 2/49-51 73 75 188 190 3/193-195 288-289 "
     "291 756-758 764 766 4/770-771 1160-1161 1189 1191 1199 3036-3038 3070" },
   { 1, 184.662, 37.615, 150.0,
     "0/0-3 5-7 9-10 1/17-19 32-33 35 44 46-47" },
   { 4, 51.269, 50.229, 170.0,
     "0/0-9 11 1/43 2/160 162 165 167-168 170-171 3/644-646 656-657 659 "
     "667 676 678-679 4/2588-2590 2608-2610 2616 2618 2632-2633 2661 2663 "
     "2667 2709-2711" },
   { 1, 242.749, 25.407, 150.0,
     "0/0-4 6-7 9-11 1/21 23 32 34" },
   { 4, 174.805, 1.291, 170.0,
     "0/0-3 5-11 2/64-66 68-69 71-72 74-75 77-79 3/268-269 280-281 283 305 "
     "307 4/1080 1128-1129 1131 1176 1178-1179 1182-1183 1217 1225-1227" },
};

static const Poly polys[] = {
   { 3, 4, { 356.7, 17.2, 18.7, 22.9, 25.3, 3.1, 0.0, 0.0 },
     "3/283-284 286 304-305 307-310" },
   { 4, 4, { 356.7, 17.2, 18.7, 22.9, 25.3, 3.1, 0.0, 0.0 },
     "3/283 286 305 308-309 4/1127 1131 1138-1139 1150 1216-1217 1219 "
     "1228-1229 1231 1240-1241" },
   { 3, 3, { 200.7, 55.1, 100.5, 50.2, 0.0, 90.0 },
     "2/30-31 3/117-119 175 186-187" },
   { 3, 3, { 190.2, -48.8, 300.1, -40.3, 90.0, -89.0 },
     "2/160-162 3/652 656 706 712 714-715 736-738" },
   { 4, 3, { 180.0, 0.0, 201.3, -14.2, 170.9, -25.6 },
     "3/391 396-397 4/1551 1561-1563 1569 1572-1573 1575 1593 1596-1597 "
     "1599 1602 1608" },
   { 4, 4, { 350.0, 10.0, 10.0, 10.0, 10.0, -10.0, 350.0, -10.0 },
     "3/271 282 293 304 4/1077-1079 1081-1083 1122 1134 1169 1181 "
     "1220-1222 1224-1226" },
   { 3, 4, { 297.267, 60.0, 207.267, 60.0, 117.267, 60.0, 27.267, 60.0 },
     "2/14-15 30-31 46-47 62-63 3/55 119 183 247" },
   { 3, 3, { 140.508, -60.0, 260.508, -60.0, 20.508, -60.0 },
     "2/128 144 3/520 522 580 582 640-642 644-645 656 704-706 712" },
   { 3, 3, { 119.168, -7.346, 67.857, -20.774, 81.401, 28.852 },
     "2/83 86 92 3/330-331 338-339 353 356-357 359 376-377" },
   { 3, 5, { 279.103, 15.996, 243.624, 23.955, 225.007, 55.249, 299.56,
       71.453, 309.149, 37.631 },
     "2/36-37 39 58-59 127 3/153 155 177 180-181 226-228 230-231 248 250 "
     "502-503 507" },
   { 3, 3, { 109.978, -64.186, 319.575, -62.602, 34.585, -36.609 },
     "2/129-131 3/513-515 544-545 548 586 708-710" },
   { 3, 5, { 167.87, 24.255, 196.987, 49.904, 231.74, 28.12, 218.74,
       -3.689, 184.501, -5.84 },
     "2/34 40 102-103 109 111 3/130-131 140 142 164 166 168-170 172 403 "
     "405-407 432-433 435 441 443" },
};

static const char *cds_ascii =
   "2/145-147 152-153 3/516-517 528-529 577-579 592 594-595 600-602 "
   "616-618 620 624";
static const char *cds_json =
   "{ \"2\": [145, 146, 147, 152, 153], \"3\": [516, 517, 528, 529, 577, "
   "578, 579, 592, 594, 595, 600, 601, 602, 616, 617, 618, 620, 624] }";

static int status = 0;
static const double d2r = 0.017453292519943295;

static void fail( const char *text ) {
   if( status == 0 ) status = 1;
   printf( "%s\n", text );
}

/* Make a Moc at the given order from a string MOC, or from the cells that
   the string MOC does not hold if "negate" is non-zero. */
static AstMoc *moc_from_string( int order, const char *text, int negate ) {
   AstMoc *moc = astMoc( "MaxOrder=%d", order );
   int json;
   astAddMocString( moc, AST__OR, negate, order, strlen( text ), text, &json );
   return moc;
}

/* Is a position, given in degrees, inside a Moc? */
static int contains( AstMoc *moc, double lon, double lat ) {
   double x = lon*d2r, y = lat*d2r, xout, yout;
   astTran2( moc, 1, &x, &y, 1, &xout, &yout );
   return astOK && xout != AST__BAD && yout != AST__BAD;
}

/* Print a Moc as a string MOC, to show what was got when a check fails. */
static void show( const char *label, AstMoc *moc ) {
   static char buf[ 8192 ];
   size_t size = 0;
   astGetMocString( moc, 0, sizeof( buf ) - 1, buf, &size );
   if( astOK ) {
      buf[ size ] = 0;
      printf( "   %s: %s\n", label, buf );
   } else {
      astClearStatus;
   }
}

/* AST must put each position in one of the cells that cdshealpix puts it
   in. */
static void check_positions( void ) {
   int n = (int)( sizeof( positions )/sizeof( positions[ 0 ] ) );
   int i, k;

   for( i = 0; i < n && status == 0; i++ ) {
      const Position *p = positions + i;
      AstMoc *in, *out;
      char text[ 128 ];
      int nc;

      astBegin;
      nc = snprintf( text, sizeof text, "%d/", p->order );
      for( k = 0; k < p->ncell; k++ ) {
         nc += snprintf( text + nc, sizeof text - (size_t) nc,
                         "%" PRId64 " ", p->cell[ k ] );
      }

      in = moc_from_string( p->order, text, 0 );
      out = moc_from_string( p->order, text, 1 );
      if( !astOK ) {
         fail( "testmoccds: could not make the Mocs for a position" );
      } else if( !contains( in, p->lon, p->lat ) ) {
         printf( "testmoccds: order %d: (%.17g,%.17g) is not in a Moc of "
                 "the cells that hold it (%s)\n", p->order, p->lon, p->lat,
                 text );
         fail( "testmoccds: a position is not in its cell" );
      } else if( contains( out, p->lon, p->lat ) ) {
         printf( "testmoccds: order %d: (%.17g,%.17g) is in a Moc of the "
                 "cells that do not hold it (all but %s)\n", p->order,
                 p->lon, p->lat, text );
         fail( "testmoccds: a position is in the wrong cell" );
      }
      astEnd;
   }
}

/* Compare the Moc of a Region with the expected string MOC. */
static void check_region( AstRegion *region, int order, const char *expected,
                          const char *what ) {
   AstMoc *moc, *want;

   moc = astMoc( "MaxOrder=%d,MinOrder=%d", order, order );
   astAddRegion( moc, AST__OR, region );
   want = moc_from_string( order, expected, 0 );
   if( !astOK ) {
      printf( "testmoccds: %s: error %d\n", what, status );
      astClearStatus;
      fail( "testmoccds: could not make the Moc of a Region" );
   } else if( !astEqual( moc, want ) ) {
      printf( "testmoccds: %s: wrong cells\n", what );
      show( "got     ", moc );
      show( "expected", want );
      fail( "testmoccds: the Moc of a Region has the wrong cells" );
   }
}

/* The Moc of each Circle must hold the cells whose centre is inside it. */
static void check_cones( void ) {
   int n = (int)( sizeof( cones )/sizeof( cones[ 0 ] ) );
   int i;

   for( i = 0; i < n && status == 0; i++ ) {
      const Cone *c = cones + i;
      char what[ 128 ];
      double centre[ 2 ], radius;

      astBegin;
      centre[ 0 ] = c->lon*d2r;
      centre[ 1 ] = c->lat*d2r;
      radius = c->radius*d2r;
      snprintf( what, sizeof what, "order %d Circle of radius %g at (%g,%g)",
                c->order, c->radius, c->lon, c->lat );
      check_region( (AstRegion *) astCircle( astSkyFrame( "System=ICRS" ), 1,
                                             centre, &radius, NULL, " " ),
                    c->order, c->moc, what );
      astEnd;
   }
}

/* The Moc of each Polygon must hold the cells whose centre is inside it. */
static void check_polys( void ) {
   int n = (int)( sizeof( polys )/sizeof( polys[ 0 ] ) );
   int i, k;

   for( i = 0; i < n && status == 0; i++ ) {
      const Poly *p = polys + i;
      char what[ 128 ];
      double points[ 10 ];

      astBegin;
      for( k = 0; k < p->nvert; k++ ) {
         points[ k ] = p->vert[ 2*k ]*d2r;
         points[ p->nvert + k ] = p->vert[ 2*k + 1 ]*d2r;
      }
      snprintf( what, sizeof what, "order %d Polygon %d starting at (%g,%g)",
                p->order, i, p->vert[ 0 ], p->vert[ 1 ] );
      check_region( (AstRegion *) astPolygon( astSkyFrame( "System=ICRS" ),
                                              p->nvert, p->nvert, points,
                                              NULL, " " ),
                    p->order, p->moc, what );
      astEnd;
   }
}

/* The string and JSON MOCs written by cds-moc for the cells at order 3
   whose centre is within 25 degrees of (123.4,-56.7) must be read as the
   Moc of that Circle.  The string MOC that AST writes for it must separate
   its cells with spaces, as the MOC recommendation requires and as cds-moc
   expects, and must be read back as the same Moc. */
static void check_text( void ) {
   static char buf[ 2048 ];
   AstMoc *moc, *back, *from_ascii, *from_json;
   double centre[ 2 ] = { 123.4*0.017453292519943295,
                          -56.7*0.017453292519943295 };
   double radius = 25.0*0.017453292519943295;
   size_t size = 0;
   int json;

   if( status != 0 ) return;
   astBegin;

   moc = astMoc( "MaxOrder=3,MinOrder=3" );
   astAddRegion( moc, AST__OR, astCircle( astSkyFrame( "System=ICRS" ), 1,
                                          centre, &radius, NULL, " " ) );

   from_ascii = astMoc( " " );
   astAddMocString( from_ascii, AST__OR, 0, -1, strlen( cds_ascii ),
                    cds_ascii, &json );
   if( astOK && ( json || !astEqual( moc, from_ascii ) ) ) {
      show( "got     ", from_ascii );
      show( "expected", moc );
      fail( "testmoccds: the string MOC written by cds-moc was misread" );
   }

   from_json = astMoc( " " );
   astAddMocString( from_json, AST__OR, 0, -1, strlen( cds_json ), cds_json,
                    &json );
   if( astOK && ( !json || !astEqual( moc, from_json ) ) ) {
      show( "got     ", from_json );
      show( "expected", moc );
      fail( "testmoccds: the JSON MOC written by cds-moc was misread" );
   }

   astGetMocString( moc, 0, sizeof( buf ) - 1, buf, &size );
   if( astOK ) {
      buf[ size ] = 0;
      if( strchr( buf, ',' ) ) {
         printf( "testmoccds: AST wrote the string MOC '%s'\n", buf );
         fail( "testmoccds: a string MOC has a comma between cells" );
      }
      back = astMoc( " " );
      astAddMocString( back, AST__OR, 0, -1, size, buf, &json );
      if( astOK && !astEqual( moc, back ) ) {
         fail( "testmoccds: a string MOC was not read back as written" );
      }
   }

   if( !astOK ) {
      astClearStatus;
      fail( "testmoccds: error checking the serialisations of a MOC" );
   }
   astEnd;
}

int main( void ) {
   astWatch( &status );

   check_positions();
   check_cones();
   check_polys();
   check_text();

   if( status == 0 ) {
      printf( "All Moc checks against the CDS HEALPix library passed\n" );
   } else {
      printf( "Moc checks against the CDS HEALPix library failed\n" );
   }
   return status != 0;
}
