#line 3367 "/Users/goshng/Dropbox/Documents/Projects/Peach/tux/noweb/main.nw"
// Example code for EZLOGGER macros
#include "ezlogger_headers.hpp"

void some_funcxx()
{
        EZLOGGERFUNCTRACKER;
}

void some_funcx()
{
        EZLOGGERFUNCTRACKER;
        some_funcxx();
}

void some_func6()
{
        EZLOGGERFUNCTRACKER;
        EZLOGGERDISPLAY_STACK;
}

void some_func5(int &x)
{
        EZLOGGERFUNCTRACKER;
        --x;
        if (x > 0)
                some_func5(x); //test recursion
        else
                some_func6();
}

void some_func4()
{
        EZLOGGERFUNCTRACKER;
        int x = 3;
        some_func5(x);
}

void some_func3()
{
        EZLOGGERFUNCTRACKER;
        some_funcx();
        some_funcxx();
        some_func4();
}

void some_func2()
{
        EZLOGGERFUNCTRACKER;
        some_funcx();
        some_funcxx();
        some_func3();
        some_funcxx();
}

void some_func1()
{
        EZLOGGERFUNCTRACKER;
        some_funcx();
        some_func2();
        some_funcxx();
}

int main(int argc, char**argv)
{
        axter::ezlogger<>::set_verbosity_level_tolerance(axter::log_very_rarely);
        EZLOGGERFUNCTRACKER;
        int ReturnValue = 99;
        EZLOGGER_PRG_MAIN_ARG(argc, argv);
        EZDBGONLYLOGGER_PRG_MAIN_ARG(argc, argv);
        EZLOGGERVL_PRG_MAIN_ARG(axter::log_often, argc, argv);
        int i = 123;
        std::string somedata = "Hello World";
        EZLOGGER(i);
        EZDBGONLYLOGGER(i);
        EZLOGGERVL(axter::log_often)(i);

        EZLOGGERVAR(somedata);
        EZDBGONLYLOGGERVAR(somedata);
        EZLOGGERVLVAR(axter::log_often, somedata);
        
        bool SomeConditionVar = true;
        EZLOGGERVAR(SomeConditionVar == false);
        EZDBGONLYLOGGERVAR(SomeConditionVar == false);
        EZLOGGERVLVAR(axter::log_often, SomeConditionVar == true);

        EZLOGGERVLVARIFY(axter::log_often, SomeConditionVar == false);

        EZLOGGERSTREAM << somedata << " " << i << std::endl;
        EZLOGGERSTREAM << somedata << " next line " << i << std::endl;
        EZLOGGERSTREAM2(std::cerr) << somedata << " next line " << i << std::endl;
        EZDBGONLYLOGGERSTREAM << somedata << " " << i << std::endl;
        EZDBGONLYLOGGERSTREAM << somedata << " next line " << i << std::endl;
        EZLOGGERVLSTREAM(axter::log_often) << somedata << " " << i << std::endl;
        // EZLOGGERVLSTREAM(axter::levels(axter::log_often, axter::warn, __FUNCSIG__ /*or GNU PRETTY_FUNCTION*/,"Xyz Facility")) << somedata << " " << i << std::endl;

        EZLOGGERPRINT("i = %i and somedata = %s", i, somedata.c_str());
        EZDBGONLYLOGGERPRINT("i = %i and somedata = %s", i, somedata.c_str());
        EZLOGGERVLPRINT(axter::log_often)("i = %i and somedata = %s", i, somedata.c_str());
        //Alternative method
        EZLOGGERVL(axter::log_often).cprint("i = %i and somedata = %s", i, somedata.c_str());
        EZLOGGER.cprint("i = %i and somedata = %s", i, somedata.c_str());

        if (1)
        {
                EZLOGGERMARKER;
                EZDBGONLYLOGGERMARKER;
                EZLOGGERVLMARKER(axter::log_often);
        }

        some_func1();

        return EZLOGGERVAR(ReturnValue);
}
