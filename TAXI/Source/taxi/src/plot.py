#!/usr/bin/python

from optparse import OptionParser
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import taxi_pickups
from const import Const

class Plotter(object):
    '''
    The Plotter class is for making plots describing the entire dataset.
    Usage:
        python plot.py --local
    '''

    def __init__(self, database, table_name):
        self.db = database
        query_string = ('SELECT * FROM ' + table_name)
        self.data = self.db.execute_query(query_string)
        matplotlib.rcParams.update({'font.size': 18})

    def plotNumPickups(self):
        '''
        Plot a histogram showing the distribution of true numbers of pickups.
        '''
        # Get data into array.
        num_pickups = []
        for row in self.data:
            num_pickups.append(row['num_pickups'])

        # Plot histogram.
        num_bins = 100
        plt.hist(num_pickups, num_bins, alpha=0.5)
        plt.yscale('log')
        plt.legend()

        # Label histogram.
        plt.title('Histogram of the number of taxi pickups')
        plt.xlabel('Number of taxi pickups in any zone and hour-long ' +
                   'time slot')
        plt.ylabel('Frequency')

        plt.show()

    def plotNumPickupsByDay(self):
        '''
        Plot a histogram showing the distribution of true number of pickups.
        by day of week. NOTE: This plot may be misleading, since there are
        more Tuesdays/Wednesdays/Thursdays in this month (Jan 2013) than
        other days.
        '''
        # Get data into arrays.
        num_days_in_week = 7
        days_in_week = np.array(range(num_days_in_week))
        num_pickups = [0]*num_days_in_week

        for row in self.data:
            day = row['start_datetime'].weekday()
            num_pickups[day] += row['num_pickups']

        # Plot bar chart.
        width = 0.7
        fig, ax = plt.subplots()
        ax.bar(days_in_week + width / 2.0, num_pickups, width=width, alpha=0.5)
        ax.set_xticks(days_in_week + width)

        # Label bar chart.
        plt.title('Number of taxi pickups by day of week')
        plt.xlabel('Day of week')
        plt.ylabel('Number of taxi pickups')
        ax.set_xticklabels(('Mon', 'Tues', 'Wed', 'Thurs', 'Fri', 'Sat', 'Sun'))

        plt.show()


    def plotNumPickupsByHour(self):
        '''
        Plot a histogram showing the distribution of true number of pickups.
        '''
        # Get data into arrays.
        num_hours = 24
        hours_in_day = np.array(range(num_hours))
        num_pickups = [0]*num_hours

        for row in self.data:
            hour = row['start_datetime'].hour
            num_pickups[hour] += row['num_pickups']

        # Plot bar chart.
        width = 1
        fig, ax = plt.subplots()
        ax.bar(hours_in_day + width / 2.0, num_pickups, width=width, alpha=0.5)
        ax.set_xticks(hours_in_day+width)

        # Label bar chart.
        plt.title('Number of taxi pickups by hour of day')
        plt.xlabel('Hour of day')
        plt.ylabel('Number of taxi pickups')
        x_labels = hours_in_day.copy().tolist()
        for idx, hour in enumerate(x_labels):
            suffix = 'am' if hour < 12 else 'pm'
            hour = hour % 12
            x_labels[idx] = '' if hour % 2 == 0 else str(hour) + suffix
        ax.set_xticklabels(x_labels)

        plt.show()

    def plotNumPickupsByZone(self):
        '''
        Plot a histogram showing the distribution of true number of pickups by
        zone.
        '''
        # Gather data.
        num_pickups_by_zone = {}
        num_pickups_list = []

        for row in self.data:
            zone_id = str(row['zone_id'])
            num_pickups_by_zone[zone_id] = \
                num_pickups_by_zone.get(zone_id, 0) + row['num_pickups']

        for num_pickups in num_pickups_by_zone.values():
            num_pickups_list.append(num_pickups)

        # Plot histogram.
        num_bins = 20
        plt.hist(num_pickups_list, num_bins, alpha=0.5)
        plt.yscale('log')

        # Label histogram.
        plt.title('Histogram of the number of taxi pickups by zone')
        plt.xlabel('Total number of taxi pickups during four months')
        plt.ylabel('Number of zones')

        plt.show()

def main(args):
    parser = OptionParser()
    parser.add_option('-l', '--local',
                      action='store_true', dest='local', default=False,
                      help='use local as opposed to remote MySQL server')
    options, args = parser.parse_args()

    database = taxi_pickups.Database(is_local=options.local)
    plotter = Plotter(database, Const.AGGREGATED_PICKUPS)
    plotter.plotNumPickups()
    plotter.plotNumPickupsByDay()
    plotter.plotNumPickupsByHour()
    plotter.plotNumPickupsByZone()

if __name__ == '__main__':
    main(sys.argv)
