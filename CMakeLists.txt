add_executable(pico_drone
        pico_drone.c
        )

# Pull in our (to be renamed) simple get you started dependencies
target_link_libraries(pico_drone pico_stdlib hardware_i2c)

# enable usb output, disable uart output
pico_enable_stdio_usb(pico_drone 1)
pico_enable_stdio_uart(pico_drone 0)

# create map/bin/hex file etc.
pico_add_extra_outputs(pico_drone)

# add url via pico_set_program_url
example_auto_set_url(pico_drone)
