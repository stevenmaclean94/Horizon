﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MissionElements;
using HSFUniverse;
using Utilities;

namespace HSFSystem
{
    public interface ISubsystem
    {
        bool CanPerform(Event proposedEvent, Universe environment);
        bool CanExtend(Event proposedEvent, Universe environment, double evalToTime); 

        void CollectDependencyFuncs(Dependency Deps, List<string> FuncNames);

    //    void DependencyCollector();
    }
}
