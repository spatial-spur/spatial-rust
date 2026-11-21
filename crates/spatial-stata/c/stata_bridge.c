/* Stata initializes the SPI callback table through stplugin.c. */
#include "stplugin.h"

extern int rs_stata_call(int argc, char *argv[]);

STDLL stata_call(int argc, char *argv[])
{
    if (!_stata_) return 198;
    return rs_stata_call(argc, argv);
}

int st_error(const char *message)
{
    return SF_error((char *)message);
}

int st_display(const char *message)
{
    return SF_display((char *)message);
}

ST_int st_in1(void)
{
    return SF_in1();
}

ST_int st_in2(void)
{
    return SF_in2();
}

ST_int st_nobs(void)
{
    return SF_nobs();
}

ST_int st_nvar(void)
{
    return SF_nvar();
}

ST_int st_nvars(void)
{
    return SF_nvars();
}

ST_boolean st_ifobs(ST_int obs)
{
    return SF_ifobs(obs);
}

ST_boolean st_is_missing(ST_double z)
{
    return SF_is_missing(z);
}

ST_double st_missval(void)
{
    return SV_missval;
}

ST_retcode st_vdata(ST_int var, ST_int obs, ST_double *z)
{
    return SF_vdata(var, obs, z);
}

ST_retcode st_vstore(ST_int var, ST_int obs, ST_double z)
{
    return SF_vstore(var, obs, z);
}

ST_int st_row(const char *mat)
{
    return SF_row((char *)mat);
}

ST_int st_col(const char *mat)
{
    return SF_col((char *)mat);
}

ST_retcode st_mat_el(const char *mat, ST_int row, ST_int col, ST_double *z)
{
    return SF_mat_el((char *)mat, row, col, z);
}

ST_retcode st_mat_store(const char *mat, ST_int row, ST_int col, ST_double z)
{
    return SF_mat_store((char *)mat, row, col, z);
}

ST_retcode st_scal_save(const char *scal, ST_double z)
{
    return SF_scal_save((char *)scal, z);
}
