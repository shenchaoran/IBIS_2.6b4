# 封装
将程序中所有写死的输入文件路径提取到参数列表中，运行时需要配置的参数有：

|类别|列表|何时需要|描述|
| - | - | - | - |
|配置输入| ini_infile_<br>diag_infile_|每次||
|站点输入|soita_sand_<br>soita_clay_<br>vegtype_<br>surta_（陆面掩膜）<br>topo_|每次||
|气象数据|cld_mon_<br>deltat_mon_<br>prec_mon_<br>rh_mon_<br>temp_mon_<br>trange_mon_<br>wetd_mon_<br>wspd_mon_|每次|月步长|
|月异常输入|temp_danom_<br>trange_danom_<br>prec_danom_<br>cld_danom_<br>rh_danom_<br>wspd_danom_<br>wetd_danom_|当ini_infile中指定需要异常数据时||
|日均值输入|prec_fanom_<br>temp_danom_<br>trange_fanomc_<br>cld_danom_<br>wspd_fanomc_<br>sphum_fanom_<br>//prec_daily_<br>//temp_daily_<br>//trange_daily_<br>//cld_daily_<br>//wspd_daily_<br>//sphum_daily_|当ini_infile中指定需要日均值数据时||
|参数|params_can_<br>params_hyd_<br>params_soi_<br>params_veg_<br>|每次|
|年输出|out_yearly_aet_<br>out_yearly_biomass_<br>out_yearly_co2fluxes_<br>out_yearly_csoi_<br>out_yearly_disturbf_<br>out_yearly_exist_<br>out_yearly_fcover_<br>out_yearly_npp_<br>out_yearly_nsoi_<br>out_yearly_plai_<br>out_yearly_runoff_<br>out_yearly_sens_<br>out_yearly_tsoi_<br>out_yearly_vegtype0_<br>out_yearly_wsoi_<br>out_yearly_zcanopy_<br>out_yearly_sapfrac_<br>out_yearly_dummyv_<br>out_yearly_solar_<br>out_yearly_albedo_<br>out_yearly_latent_<br>out_yearly_totfall_<br>out_yearly_clitw_<br>out_yearly_csoislo_<br>out_yearly_csoipas_|指定需要数据的项目路径，不指定的不予输出|
|月输出|out_monthly_aet_<br>out_monthly_cloud_<br>out_monthly_co2ratio_<br>out_monthly_ir_<br>out_monthly_lai_<br>out_monthly_latent_<br>out_monthly_npptot_<br>out_monthly_qa_<br>out_monthly_rain_<br>out_monthly_rh_<br>out_monthly_runoff_<br>out_monthly_sens_<br>out_monthly_snod_<br>out_monthly_snof_<br>out_monthly_snow_<br>out_monthly_solar_<br>out_monthly_temp_<br>out_monthly_tsoi_<br>out_monthly_wsoi_<br>out_monthly_albedo_<br>out_monthly_dummyv_<br>|指定需要数据的项目路径，不指定的不予输出|
|日输出|out_daily_rain_<br>out_daily_cloud_<br>out_daily_rh_<br>out_daily_snow_<br>out_daily_aet_<br>out_daily_trunoff_<br>out_daily_srunoff_<br>out_daily_drainage_<br>out_daily_wsoi_<br>out_daily_wisoi_<br>out_daily_snod_<br>out_daily_snof_<br>out_daily_co2ratio_<br>out_daily_co2mic_<br>out_daily_templ_<br>out_daily_zcanopy_<br>out_daily_laicanopy_<br>|指定需要数据的项目路径，不指定的不予输出|
|诊断输出|out_diag_0_<br>out_diag_1_<br>out_diag_2_<br>out_diag_3_<br>out_diag_4_<br>out_diag_5_<br>out_diag_6_<br>out_diag_7_<br>out_diag_8_<br>out_diag_9_<br>|指定需要数据的项目路径，不指定的不予输出|
|其他输出|out_global_<br>out_vegtype_<br>out_yearsrun_<br>|指定需要数据的项目路径，不指定的不予输出|

# 编译
## makefile
- make ibis
- make clean

## 编译选项
- F77_OPTIONS：编译器，f77, gfortron, ifort...
- INCLUDE_DIRS：头文件目录
- LD_OPTIONS_NETCDF：lib目录，使用的 libnetcdff 是动态链接库(so)

# 数据准备
## 全球文件
标准全球模式需要13个输入文件，这些文件并不是原生的 netcdf 格式。我们准备了一份数据，如果你同意在应用时添加引用，联系我们获取这些数据。

## NetCDF
## 分辨率
我们用的原始数据分辨率是 0.5*0.5. 我们将其重采样为 1.0, 2,0 和 4.0 的。你可以创建你自己的空间分辨率文件。

## 13个输入文件
8个气象数据，5个站点数据

| Name | Filename | Description | Units | Source | Type |
| ---- | -------- | ----------- | ----- | ------ |-|
| cld     | cld.mon.nc    | monthly mean cloudiness                                                     | %              | CRU                 | Met |
| deltat  | deltat.mon.nc | minimum temp ever recorded at that location minus avg temp of coldest month | C              | Oregon              | Met |
| prec    | prec.mon.nc   | monthly mean precipitation rate                                             | mm/day         | CRU                 | Met |
| rh      | rh.mon.nc     | monthly mean relative humidity                                              | %              | CRU                 | Met |
| temp    | temp.mon.nc   | monthly mean temperature                                                    | C              | CRU                 | Met |
| trange  | trange.mon.nc | monthly mean temperature range                                              | C              | CRU                 | Met |
| wetd    | wetd.mon.nc   | mean "wet" days per month                                                   | days           | CRU                 | Met |
| wspd    | wspd.mon.nc   | monthly mean wind speed at sig=0.995                                        | m/s            | CRU                 | Met |
| sand    | soita.sand.nc | percentage of sand                                                          | %              | IGBP/CONUS          | Site |
| clay    | soita.clay.nc | percentage of clay                                                          | %              | IGBP/CONUS          | Site |
| vegtype | vegtype.nc    | initial vegetation types                                                    | indexd         | SAGE                | Site |
| surta   | surta.nc      | land mask                                                                   | 1=land,0=ocean | adapted from ETOPO5 | Site |
| topo    | topo.nc       | topography                                                                  | m              | ETOPO5              | Site |

### CRU：全球气候数据集
### OREGON：
### IGBP：全球土壤数据集
### CONUS：美国土壤数据集
### ETOP05：surta 和 topo
### SAGA：植被类型数据
### 创建自己的数据集

## 测试数据
我们提供了一些没有实际意义的测试数据，空间分布率为 10*10，用来测试编译结果。

## 真实数据
我们提供了0.5，1.0，2.0，4.0的空间分辨率的真实数据，联系我们获取。

# 运行
1. 创建输入输出文件路径
2. 将输入文件放在输入路径下。输入数据由 ies-io.f 读取
3. 编辑 compar.h 文件，里面设定了输入数据的元信息。这一部分待封装
    - nlon：
    - nlat：
    - npoi：
    - xres：
    - yres：
4. 编辑 ibis.infile。配置运行参数。
    - irestart：0（new run）；1（restart）
    - iyear0：模拟的第一年年份。重启时不能改变。
    - nrun：模拟总年数。重启时必须变成剩余的年数。比如，初始的模拟总年数是100，在第80年重启了，那么需要将该值设为21。
    - nanom：开始读取**月异常**的第一年年份。如果不想读取异常数据，将此值设置成非常大（大于 iyear0 + nrun - 1）。
    - ndprecy：开始读取**日均值**的第一年年份。如果不想读取日均值，将此值设置成非常大（大于 iyear0 + nrun -1）。
    - soilcspin：1（开启土地 spinup 加速）；0（不开启）
    - iyearout：1（输出年数据）；0（不输出）
    - imonthout：1（输出月数据）；0（不输出）
    - idailyout：1（输出日数据）；0（不输出）
    - isimveg：0（静态植被）；1（动态植被）；2（从0开始的动态植被，即初始没有植被）。注意，该值为0或1时需要输入 vegtype.nc。
    - isimfire：0（固定着火状态）；1（动态着火状态）。推荐值0.
    - isimco2：0（固定co2浓度）；1（根据多项式算出来的co2浓度）
    - co2init：初始co2浓度（mol/mol）
    - o2init：初始o2浓度（mol/mol）
    - dtime：模拟的时间步长（秒）。必须是 86400（24小时）的偶除数
    - idiag：0（不输出诊断）；1-10（输出诊断数据，并设置诊断文件的数量）
    - snorth：子区域北天门
    - ssouth：子区域南天门
    - swest：子区域西天门
    - seast：子区域东天门
5. 编辑 diag.infile，用于每隔一定的时间步长打印选择的变量。
6. 运行 ibis，在命令行中输入 `ibis` 就行了
7. 重启运行时不需要重新编译，第3条涉及的项目改变时需要重新编译。

# 参数
## wpudmax
## nsoilay

# NetCDF 文件的读写
## 两个矢量：istart, icount
在 Fortran 中，多维数组中第一维的读写最快，与在C语言中相反。

如果使用 ncdump（c编写的读写库） 读一个 NetCDF 文件，其中 mydata 变量含有4维：（long, lat, level, time）。变量...

NetCDF 的读写通过两个矢量实现：istart, icount。istart表示沿着每一维度读写的起始点：`(startIndexLong, startIndexLat, startIndexLevel, startIndexTime)`，icount表示每一维数据读取多少：`(countLong, countLat,countLevel, countTime)`。

1. 读取整个区域的一个时间步长（6th）的一个变量（2nd）：`istart = (1, 1, 2, 6); icount = (nlons, nlats, 1, 1);`
2. 读取整个区域的的第6个时间步长的9个变量：`istart = (1, 1, 1, 6); icount = (nlons, nlats, 9, 1);`
3. 三维变量中单独一点的所有时间步长内的值：`istart(ilon, ilat,1); icount = (ilon, ilat, 100)`。如果istart 和icount 被声明为4维，而它用来指向一个3维变量时，第4维将被忽略。
4. 网格的子区域中，12个时间步长的一个变量指标：`istart = (ilon, ilat, 18, itime); icount = (20, 15, 1, 12)`

## 如何读取文件
在 ies-io.f 文件的子程序 readvar 中。
`call readvar(filen,varname,name3d,istart,icount,values, alons,alats,vals3d,times,ierror)`

**输入：**
- filen：文件名
- varname：要读的变量名
- name3d：第3维非时间维变量名，如果varname变量是3维的就忽略此值。
- istart
- icount

**输出：**
- values：
- alons：
- alats：
- vals3d：
- times：
- ierror：

## 如何写文件