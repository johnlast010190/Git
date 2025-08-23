
# User settings file may be a symlink, in which case install a copy
#get_filename_component(real_settings_file ${HELYX_SETTINGS_FILE} REALPATH)
configure_file(@real_settings_file@
    @cpack_temporary_source_dir@/etc/userSettings.cmake
    COPYONLY
    USE_SOURCE_PERMISSIONS
    )
