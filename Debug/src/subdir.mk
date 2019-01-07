################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Mcoption.cpp \
../src/Sdepaths.cpp \
../src/particle.cpp 

OBJS += \
./src/Mcoption.o \
./src/Sdepaths.o \
./src/particle.o 

CPP_DEPS += \
./src/Mcoption.d \
./src/Sdepaths.d \
./src/particle.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/oliv/Documents/GitHub/topov2" -O0 -g3 -Wall -c -fmessage-length=0 -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


