#ifndef COMMON_H
#define COMMON_H


#define     VSP             vtkSmartPointer

enum OPERATE
{
    DATA_LOAD = 0,
    TOOTH_DESIGN,
    POSITIONER_DESIGN,
    PLAN_GENERATE,
    OPERATE_COUNT
};

enum ACTION_ID
{
    AC_INVALID = -1,
    AC_TO_ORDER_LIST,
    AC_LOAD_UPPER,
    AC_LOAD_LOWER,

    AC_COUNT
};


const int TOOTH_COUNT = 32;

enum TOOTH_INDEX
{
    INDEX_INVALID = 0,
    INDEX_01 = 1,
    INDEX_02,
    INDEX_03,
    INDEX_04,
    INDEX_05,
    INDEX_06,
    INDEX_07,
    INDEX_08,
    INDEX_09,
    INDEX_10,
    INDEX_11,
    INDEX_12,
    INDEX_13,
    INDEX_14,
    INDEX_15,
    INDEX_16,
    INDEX_17,
    INDEX_18,
    INDEX_19,
    INDEX_20,
    INDEX_21,
    INDEX_22,
    INDEX_23,
    INDEX_24,
    INDEX_25,
    INDEX_26,
    INDEX_27,
    INDEX_28,
    INDEX_29,
    INDEX_30,
    INDEX_31,
    INDEX_32
};

enum TOOTH_ID
{
    ID_INVALID = 0,
    ID_11 = 11,
    ID_12,
    ID_13,
    ID_14,
    ID_15,
    ID_16,
    ID_17,
    ID_18,
    ID_21 = 21,
    ID_22,
    ID_23,
    ID_24,
    ID_25,
    ID_26,
    ID_27,
    ID_28,
    ID_31 = 31,
    ID_32,
    ID_33,
    ID_34,
    ID_35,
    ID_36,
    ID_37,
    ID_38,
    ID_41 = 41,
    ID_42,
    ID_43,
    ID_44,
    ID_45,
    ID_46,
    ID_47,
    ID_48
};

extern TOOTH_INDEX ID_TO_INDEX[49];
extern TOOTH_ID INDEX_TO_ID[33];

#endif // COMMON_H
