################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../assign3.o \
../intersection_point.o \
../ray.o \
../vector.o 

CPP_SRCS += \
../assign3.cpp \
../intersection_point.cpp \
../ray.cpp \
../vector.cpp 

OBJS += \
./assign3.o \
./intersection_point.o \
./ray.o \
./vector.o 

CPP_DEPS += \
./assign3.d \
./intersection_point.d \
./ray.d \
./vector.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


