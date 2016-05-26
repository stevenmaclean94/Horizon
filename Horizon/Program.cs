﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Xml;
using System.Text;
using HSFScheduler;
using Utilities;
using MissionElements;
using UserModel;
using HSFUniverse;
using HSFSubsystem;
using HSFSystem;

namespace Horizon
{
    public class Program
    {
        static void Main(string[] args)
        {
            // Get the input filenames
            //string simulationInputFilePath = args[1];
            //string targetDeckFilePath = args[2];
            //string modelInputFileName = args[3];
            //string outputPath = args[4];
            var simulationInputFilePath = @"..\..\..\SimulationInput.XML"; // @"C:\Users\admin\Documents\Visual Studio 2015\Projects\Horizon-Simulation-Framework\Horizon_v2_3\io\SimulationInput.XML";
            var targetDeckFilePath = @"..\..\..\v2.2-300targets.xml";
            var modelInputFilePath = @"..\..\..\DSAC_Static.xml";

            var outputFileName = string.Format("output-{0:yyyy-MM-dd-hh-mm-ss}-*", DateTime.Now);
            var outputPath = @"..\..\..\";
            var txt = ".txt";
            string[] fileNames = System.IO.Directory.GetFiles(outputPath, outputFileName, System.IO.SearchOption.TopDirectoryOnly);
            double number = 0;
            foreach (var fileName in fileNames)
            {
                char version = fileName[fileName.Length - txt.Length-1];
                if(number < Char.GetNumericValue(version))
                    number = Char.GetNumericValue(version);
            }
            number++;
            outputFileName = outputFileName.Remove(outputFileName.Length - 1) + number;
            outputPath += outputFileName + txt;
            // Find the main input node from the XML input files
            var XmlDoc = new XmlDocument();
            XmlDoc.Load(simulationInputFilePath);
            XmlNodeList simulationInputXMLNodeList = XmlDoc.GetElementsByTagName("SCENARIO");
            var XmlEnum = simulationInputXMLNodeList.GetEnumerator();
            XmlEnum.MoveNext();
            var simulationInputXMLNode = (XmlNode)XmlEnum.Current;
            var scenarioName = simulationInputXMLNode.Attributes["scenarioName"].InnerXml;
            Console.Write("EXECUITING SCENARIO: ");
            Console.WriteLine(scenarioName);

            // Load the simulation parameters from the XML simulation input file
            XmlNode simParametersXMLNode = simulationInputXMLNode["SIMULATION_PARAMETERS"];
            bool simParamsLoaded = SimParameters.LoadSimParameters(simParametersXMLNode, scenarioName);

            // Load the scheduler parameters defined in the XML simulation input file
            XmlNode schedParametersXMLNode = simulationInputXMLNode["SCHEDULER_PARAMETERS"];

            bool paramsLoaded = SchedParameters.LoadSchedParameters(schedParametersXMLNode);
            Scheduler systemScheduler = new Scheduler();
            //MultiThreadedScheduler systemScheduler = new MultiThreadedScheduler();

            // Load the target deck into the targets list from the XML target deck input file
            XmlDoc.Load(targetDeckFilePath);
            XmlNodeList targetDeckXMLNodeList = XmlDoc.GetElementsByTagName("TARGETDECK");
            int numTargets = XmlDoc.GetElementsByTagName("TARGET").Count;
            XmlEnum = targetDeckXMLNodeList.GetEnumerator();
            XmlEnum.MoveNext();
            var targetDeckXMLNode = (XmlNode)XmlEnum.Current;
            Stack<Task> systemTasks = new Stack<Task>();
            bool targetsLoaded = Task.loadTargetsIntoTaskList(targetDeckXMLNode, systemTasks);
            Console.WriteLine("Initial states set");

            // Find the main model node from the XML model input file
            XmlDoc.Load(modelInputFilePath);
            XmlNodeList modelXMLNodeList = XmlDoc.GetElementsByTagName("MODEL");
            XmlEnum = modelXMLNodeList.GetEnumerator();
            XmlEnum.MoveNext();
            var modelInputXMLNode = (XmlNode)XmlEnum.Current;

            // Load the environment. First check if there is an ENVIRONMENT XMLNode in the input file
            Universe SystemUniverse = null;
            foreach (XmlNode node in modelInputXMLNode.ChildNodes)
            {
                if (node.Attributes["ENVIRONMENT"] != null)
                {
                    // Create the Environment based on the XMLNode
                    SystemUniverse = new Universe(node);
                }
            }
            if (SystemUniverse == null)
                SystemUniverse = new Universe();

            //Create singleton dependency dictionary
            Dependency dependencies = Dependency.Instance;

            // Initialize List to hold assets and subsystem nodes
            List<Asset> assetList = new List<Asset>();
            List<Subsystem> subList = new List<Subsystem>();

            // Maps used to set up preceeding nodes
            Dictionary<ISubsystem, XmlNode> subsystemXMLNodeMap = new Dictionary<ISubsystem, XmlNode>();
            Dictionary<string, Subsystem> subsystemMap = new Dictionary<string, Subsystem>();
            List<KeyValuePair<string, string>> dependencyMap = new List<KeyValuePair<string, string>>();
            List<KeyValuePair<string, string>> dependencyFcnMap = new List<KeyValuePair<string, string>>();
            // Dictionary<string, ScriptedSubsystem> scriptedSubNames = new Dictionary<string, ScriptedSubsystem>();

            // Create Constraint list 
            List<Constraint> constraintsList = new List<Constraint>();

            //Create Lists to hold all the initial condition and dependency nodes to be parsed later
            List<XmlNode> ICNodes = new List<XmlNode>();
            List<XmlNode> DepNodes = new List<XmlNode>();
            SystemState initialSysState = new SystemState();

            // Set up Subsystem Nodes, first loop through the assets in the XML model input file
            foreach (XmlNode childNodeAsset in modelInputXMLNode.ChildNodes)
            {
                if (childNodeAsset.Name.Equals("PYTHON"))
                {
                    // Loop through all the of the file nodes -- TODO (Morgan) What other types of things might be scripted
                    foreach (XmlNode fileXmlNode in childNodeAsset.ChildNodes)
                    {
                        // If scripting is enabled, parse the script file designated by the attribute
                            // Parse script file if the attribute exists
                            if (fileXmlNode.ChildNodes[0].Name.Equals("EOMS_FILE"))
                            {
                                string fileName = fileXmlNode.ChildNodes[0].Attributes["src"].Value.ToString();
                                ScriptedEOMS eoms = new ScriptedEOMS(fileName);
                            }
                        }
                    }
                if (childNodeAsset.Name.Equals("ASSET"))
                {
                    Asset asset = new Asset(childNodeAsset);
                    assetList.Add(asset);
                    // Loop through all the of the ChildNodess for this Asset
                    foreach (XmlNode childNode in childNodeAsset.ChildNodes)
                    {
                        // Get the current Subsystem XML Node, and create it using the SubsystemFactory
                        if (childNode.Name.Equals("SUBSYSTEM"))
                        {  //is this how we want to do this?
                            // Check if the type of the Subsystem is scripted, networked, or other
                            string subName = SubsystemFactory.GetSubsystem(childNode, dependencies, asset, subsystemMap);
                            foreach (XmlNode ICorDepNode in childNode.ChildNodes)
                            {
                                if(ICorDepNode.Name.Equals("IC"))
                                    ICNodes.Add(ICorDepNode);
                                if (ICorDepNode.Name.Equals("DEPENDENCY"))
                                {
                                    string depSubName = "", depFunc = "";
                                    depSubName = Subsystem.parseNameFromXmlNode(ICorDepNode, asset.Name) ;
                                    dependencyMap.Add(new KeyValuePair<string, string>(subName, depSubName));

                                    if (ICorDepNode.Attributes["fcnName"] != null)
                                    {
                                        depFunc = ICorDepNode.Attributes["fcnName"].Value.ToString();
                                        dependencyFcnMap.Add(new KeyValuePair<string, string>(subName, depFunc));
                                    }
                                }  
                            }
                        }
                        //Create a new Constraint
                        if (childNode.Name.Equals("CONSTRAINT"))
                        {
                            constraintsList.Add(ConstraintFactory.getConstraint(childNode, subsystemMap, asset));
                        }
                    }
                    if (ICNodes.Count > 0)
                        initialSysState.Add(SystemState.setInitialSystemState(ICNodes, asset));
                    ICNodes.Clear();
                }
            }
            foreach (KeyValuePair<string, Subsystem> sub in subsystemMap)
            {
                subList.Add(sub.Value);
            }
            Console.WriteLine("Subsystems and Constraints Loaded");

            //Add all the dependent subsystems to the dependent subsystem list of the subsystems
            foreach (KeyValuePair<string, string> depSubPair in dependencyMap)
            {
                Subsystem subToAddDep, depSub;
                subsystemMap.TryGetValue(depSubPair.Key, out subToAddDep);
                subsystemMap.TryGetValue(depSubPair.Value, out depSub);
                subToAddDep.DependentSubsystems.Add(depSub);
            }

            //give the dependency functions to all the subsytems that need them
            foreach (KeyValuePair<string, string> depFunc in dependencyFcnMap)
            {
                Subsystem subToAddDep;
                subsystemMap.TryGetValue(depFunc.Key, out subToAddDep);
                subToAddDep.SubsystemDependencyFunctions.Add(depFunc.Value, dependencies.getDependencyFunc(depFunc.Value));
            }
            Console.WriteLine("Dependencies Loaded");

            //Need to make this parse Xml
            Evaluator schedEvaluator = new TargetValueEvaluator(dependencies);

            SystemClass simSystem = new SystemClass(assetList, subList, constraintsList, SystemUniverse);

            if (simSystem.checkForCircularDependencies())
                throw new NotFiniteNumberException("System has circular dependencies! Please correct then try again.");

            //Subsystem power;
            //Event myEvent = new Event(new Dictionary<Asset, Task>(), new SystemState());
            //subsystemMap.TryGetValue("asset1.power", out power);
            //Delegate func;

            //power.SubsystemDependencyFunctions.TryGetValue("PowerfromADCS", out func);
            //HSFProfile<double> temp = (HSFProfile<double>)func.DynamicInvoke(myEvent);
            Scheduler scheduler = new Scheduler();
            List<SystemSchedule> schedules = scheduler.GenerateSchedules(simSystem, systemTasks, initialSysState, schedEvaluator);
            // Evaluate the schedules and set their values
            foreach (SystemSchedule systemSchedule in schedules)
                systemSchedule.ScheduleValue = schedEvaluator.Evaluate(systemSchedule);

            // Sort the sysScheds by their values
            schedules.Sort((x, y) => x.ScheduleValue.CompareTo(y.ScheduleValue));
            schedules.Reverse();
            double maxSched = schedules[0].ScheduleValue;
            int i = 0;
            //Morgan's Way
            using (StreamWriter sw = File.CreateText(outputPath))
            {
                foreach (SystemSchedule sched in schedules)
                {
                    sw.WriteLine("Schedule Number: " + i + "Schedule Value: " + schedules[i].ScheduleValue);
                    foreach (var eit in sched.AllStates.Events)
                    {
                        if (eit.Tasks.Values.GetType().Equals(TaskType.COMM))
                        {
                            Console.WriteLine("Schedule {0} contains Comm task", i);
                        }
                        if (i < 5)
                        { //just compare the first 5 schedules for now
                            sw.WriteLine(eit.ToString());
                        }
                    }
                    i++;
            }
            Console.WriteLine(maxSched);
            }

            ////Mehiel's way
            string stateDataFilePath = @"..\..\..\" + string.Format("output-{0:yyyy-MM-dd-hh-mm-ss}", DateTime.Now);
            SystemSchedule.WriteSchedule(schedules[0], stateDataFilePath);

            var csv = new StringBuilder();
            csv.Clear();
            foreach (var asset in simSystem.Assets)
            {
                File.WriteAllText(@"..\..\..\" + asset.Name + "_dynamicStateData.csv", asset.AssetDynamicState.ToString());
            }

            //   Console.ReadKey();
     

                // *********************************Output selected data*************************************
             //   bool schedOutput = dataOut.writeAll(schedules, simSystem);
            // ******************************************************************************************
        }
    }
}
